#pragma once

#include <deque>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>

template < typename T >
class QueueItemInfo {
public:
  static size_t Count( const T& t ) {
    return 1;
  }
};

template < class Worker, class QueueItem, typename... Args >
class WorkerQueue {
public:
  using OnProcessedCallback =
    std::function< void( const size_t, const size_t ) >;

  WorkerQueue( const int numWorkers = 1, Args... args )
      : mStop( false ), mWorkingCount( 0 ), mTotalEnqueued( 0 ),
        mTotalProcessed( 0 ) {
    auto actualWorkers =
      numWorkers <= 0 ? std::thread::hardware_concurrency() : numWorkers;

    for( int i = 0; i < actualWorkers; i++ ) {
      mWorkers.push_back( std::thread(
        [this]( Args&&... args ) {
          this->WorkerLoop( std::forward< Args >( args )... );
        },
        args... ) );
    }
  }

  ~WorkerQueue() {
    mStop = true;
    mCondition.notify_all();
    for( auto& worker : mWorkers ) {
      if( worker.joinable() ) {
        worker.join();
      }
    }
  }

  void Enqueue( QueueItem& queueItem ) {
    {
      std::unique_lock< std::mutex > lock( mQueueMutex );
      mTotalEnqueued += QueueItemInfo< QueueItem >::Count( queueItem );
      mQueue.push( std::move( queueItem ) );
    }

    // Send one worker to work
    mCondition.notify_one();
  }

  bool Done() const {
    return mWorkingCount == 0 && mQueue.empty();
  }

  void WaitTillDone() {
    while( !Done() ) {
      std::this_thread::sleep_for( std::chrono::milliseconds( 50 ) );
    }
  }

  void OnProcessed( const OnProcessedCallback& callback ) {
    mProcessedCallbacks.push_back( callback );
  }

private:
  std::deque< std::thread > mWorkers;

  std::condition_variable mCondition;
  std::mutex              mQueueMutex;
  std::atomic< bool >     mStop;
  std::atomic< int >      mWorkingCount;

  std::queue< QueueItem > mQueue;

  size_t                            mTotalEnqueued;
  size_t                            mTotalProcessed;
  std::deque< OnProcessedCallback > mProcessedCallbacks;

  void WorkerLoop( Args&&... args ) {
    QueueItem queueItem;
    Worker    worker( std::forward< Args >( args )... );

    while( true ) {
      { // acquire lock
        std::unique_lock< std::mutex > lock( mQueueMutex );

        while( !mStop && mQueue.empty() )
          mCondition.wait( lock );

        if( mStop )
          break;

        queueItem = std::move( mQueue.front() );
        mQueue.pop();

        mWorkingCount++;
      } // release lock

      worker.Process( queueItem );

      { // acquire lock
        std::unique_lock< std::mutex > lock( mQueueMutex );
        mTotalProcessed += QueueItemInfo< QueueItem >::Count( queueItem );
        mWorkingCount--;

        for( auto& cb : mProcessedCallbacks ) {
          cb( mTotalProcessed, mTotalEnqueued );
        }
      } // release lock
    }
  }
};
