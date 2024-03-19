# List of memory cgroups that can be used for leak-y targets (subprocesses such as R
# and Julia are greedy and will not GC often enough). cgroup memory limit in
# bytes can be set to <8GB (we used 7864320000 B) for 4 cgroups each, or set to NULL to disable this
# feature.
memory_cgroups <- NULL