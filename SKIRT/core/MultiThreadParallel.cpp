/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MultiThreadParallel.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

MultiThreadParallel::MultiThreadParallel(int threadCount)
{
    // Cache the number of threads
    _threadCount = threadCount;

    // Remember the ID of the current thread
    _parentThread = std::this_thread::get_id();

    // Initialize shared data members and launch threads in a critical section
    {
        std::unique_lock<std::mutex> lock(_mutex);

        // Initialize shared data members
        _numChunks = 0;
        _chunkSize = 0;
        _maxIndex = 0;
        _active.assign(threadCount, true);
        _exception = nullptr;
        _terminate = false;
        _next = 0;

        // Create the extra parallel threads with one-based index (parent thread has index zero)
        for (int index = 1; index < threadCount; index++)
        {
            _threads.push_back(std::thread(&MultiThreadParallel::run, this, index));
        }
    }

    // Wait until all parallel threads are ready
    waitForThreads();
}

////////////////////////////////////////////////////////////////////

MultiThreadParallel::~MultiThreadParallel()
{
    // Ask the parallel threads to exit in a critical section
    {
        std::unique_lock<std::mutex> lock(_mutex);

        // Ask the parallel threads to exit
        _terminate = true;
        _conditionExtra.notify_all();
    }

    // Wait for them to do so
    for (auto& thread: _threads) thread.join();
}

////////////////////////////////////////////////////////////////////

void MultiThreadParallel::call(std::function<void(size_t,size_t)> target, size_t maxIndex)
{
    // Verify that we're being called from our parent thread
    if (std::this_thread::get_id() != _parentThread)
        throw FATALERROR("Parallel call not invoked from thread that constructed this object");

    // Intercept the pathological case where we have no work at all to avoid zero division later on
    if (!maxIndex) return;

    // Initialize shared data members and activate threads in a critical section
    {
        std::unique_lock<std::mutex> lock(_mutex);

        // Copy the target function so it can be invoked from any of the threads
        _target = target;

        // Determine the chunk size and the number of chunks
        const size_t numChunksPerThread = 8;   // empirical multiplicator to achieve acceptable load balancing
        _maxIndex = maxIndex;
        _numChunks = min(_threadCount>1 ? _threadCount*numChunksPerThread : 1, maxIndex);
        _chunkSize = maxIndex / _numChunks;
        if (_numChunks * _chunkSize < maxIndex) _numChunks++;

        // Initialize the number of active threads (i.e. not waiting for new work)
        _active.assign(_threadCount, true);

        // Clear the exception pointer
        _exception = nullptr;

        // Initialize the loop variable
        _next = 0;

        // Wake all parallel threads, if multithreading is allowed
        _conditionExtra.notify_all();
    }

    // Do some work ourselves as well
    doWork();

    // Wait until all parallel threads are done
    waitForThreads();

    // Check for and process the exception, if any
    if (_exception)
    {
        throw *_exception;  // throw by value (the memory for the heap-allocated exception is leaked)
    }
}

////////////////////////////////////////////////////////////////////

void MultiThreadParallel::run(int threadIndex)
{
    while (true)
    {
        // Wait for new work in a critical section
        {
            std::unique_lock<std::mutex> lock(_mutex);

            // Indicate that this thread is no longer doing work
            _active[threadIndex] = false;

            // Tell the main thread when all parallel threads are inactive
            if (!threadsActive()) _conditionMain.notify_all();

            // Wait for new work
            while (true)
            {
                _conditionExtra.wait(lock);

                // Check for termination request (don't bother with _active; it's no longer used)
                if (_terminate) return;

                // Check that we actually have new work
                if (_active[threadIndex]) break;
            }
        }

        // Do work as long as some is available
        doWork();
    }
}

////////////////////////////////////////////////////////////////////

void MultiThreadParallel::doWork()
{
    try
    {
        // Do work as long as some is available
        while (true)
        {
            size_t chunkIndex = _next++;             // get the next chunk index atomically
            if (chunkIndex >= _numChunks) break;     // break if no more are available

            // Invoke the target function
            _target(chunkIndex*_chunkSize, min(_chunkSize,_maxIndex-chunkIndex*_chunkSize));
        }
    }
    catch (FatalError& error)
    {
        // Make a copy of the exception
        reportException(new FatalError(error));
    }
    catch (...)
    {
        // Create a fresh exception
        reportException(new FATALERROR("Unhandled exception (not of type FatalError) in a parallel thread"));
    }
}

////////////////////////////////////////////////////////////////////

void MultiThreadParallel::reportException(FatalError* exception)
{
    // Need to lock, in case multiple threads throw simultaneously
    std::unique_lock<std::mutex> lock(_mutex);
    if (!_exception)  // only store the first exception thrown
    {
        _exception = exception;

        // Make the other threads stop by taking away their work
        _numChunks = 0;  // another thread will hopefully see either the old value, or zero
    }
}

////////////////////////////////////////////////////////////////////

void MultiThreadParallel::waitForThreads()
{
    // Wait until all parallel threads are inactive
    std::unique_lock<std::mutex> lock(_mutex);
    while (threadsActive()) _conditionMain.wait(lock);
}

////////////////////////////////////////////////////////////////////

bool MultiThreadParallel::threadsActive()
{
    // Check for active threads, skipping the parent thread with index 0
    return std::any_of(_active.begin()+1, _active.end(), [](bool flag){return flag;});
}

////////////////////////////////////////////////////////////////////