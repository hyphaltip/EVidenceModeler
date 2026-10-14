//! Thread-pool setup shared by the command-line tools.

/// Size the global rayon pool. `0` keeps rayon's default (the `RAYON_NUM_THREADS`
/// variable, else the CPUs available to the process). Call once, before any
/// parallel work; a pool that is already built is left alone.
pub fn init_thread_pool(threads: usize) {
    if threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .ok();
    }
}
