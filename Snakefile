BASE = '/oak/stanford/groups/pritch/users/jake_harold'
  # paste0('../results/mcmc/dmult_',
  #   cur_seed, '_', n_samples, '_', n_burnin, '_', thin, '.rds'))

seed_start = 42
n_chains = 4

rule all:
    input:
        expand(BASE + '/results/mcmc/dmult_{seed}_{n}_{burn}_{thin}.rds',
            seed = range(seed_start, seed_start + n_chains),
            n = 1000000, burn = 10000, thin = 20),
        expand(BASE + '/results/mcmc/dmult_flexible_bins{seed}_{n}_{burn}_{thin}.rds',
            seed = range(seed_start, seed_start + 10),
            n = 100000, burn = 50000, thin = 1),
        # expand(BASE + '/results/mcmc/dmult_{seed}_{n}_{burn}_{thin}.rds',
        #     seed = range(seed_start, seed_start + 20),
        #     n = 100000, burn = 1000, thin = 1),
        expand(BASE + '/results/mcmc/nb_{seed}_{n}_{burn}_{thin}.rds',
            seed = range(seed_start, seed_start + n_chains),
            n = 1000000, burn = 10000, thin = 20)

rule dmult_nimble_mcmc_flexible_bins:
    benchmark:
        BASE + '/benchmarks/nimble/dmult_flexible_bins{seed}_{niter}_{burn}_{thin}.txt'
    output:
        BASE + '/results/mcmc/dmult_flexible_bins{seed,\d+}_{niter}_{burn}_{thin}.rds'
    shell:
        'cd {BASE}/playground'
        ' && '
        'module load gcc/9.1.0'
        ' && '
        'module load R'
        ' && '
        'Rscript {BASE}/playground/dmult_nimble_mcmc_flexible_bins.R'
        ' {wildcards.seed}'
        ' {wildcards.niter}'
        ' {wildcards.burn}'
        ' {wildcards.thin}'

rule dmult_nimble_mcmc:
    benchmark:
        BASE + '/benchmarks/nimble/dmult_{seed}_{niter}_{burn}_{thin}.txt'
    output:
        BASE + '/results/mcmc/dmult_{seed,\d+}_{niter}_{burn}_{thin}.rds'
    shell:
        'cd {BASE}/playground'
        ' && '
        'module load gcc/9.1.0'
        ' && '
        'module load R'
        ' && '
        'Rscript {BASE}/playground/dmult_nimble_mcmc.R'
        ' {wildcards.seed}'
        ' {wildcards.niter}'
        ' {wildcards.burn}'
        ' {wildcards.thin}'


rule gp_nimble_mcmc:
    benchmark:
        BASE + '/benchmarks/nimble/gp_{seed}_{niter}_{burn}_{thin}.txt'
    output:
        BASE + '/results/mcmc/gp_{seed}_{niter}_{burn}_{thin}.rds'
    shell:
        'cd {BASE}/playground'
        ' && '
        'module load gcc/9.1.0'
        ' && '
        'module load R'
        ' && '
        'Rscript {BASE}/playground/gp_sample_mcmc.R'
        ' {wildcards.seed}'
        ' {wildcards.niter}'
        ' {wildcards.burn}'
        ' {wildcards.thin}'

rule nb_nimble_mcmc:
    benchmark:
        BASE + '/benchmarks/nimble/nb_{seed}_{niter}_{burn}_{thin}.txt'
    output:
        BASE + '/results/mcmc/nb_{seed}_{niter}_{burn}_{thin}.rds'
    shell:
        'cd {BASE}/playground'
        ' && '
        'module load gcc/9.1.0'
        ' && '
        'module load R'
        ' && '
        'Rscript {BASE}/playground/nb_sample_mcmc.R'
        ' {wildcards.seed}'
        ' {wildcards.niter}'
        ' {wildcards.burn}'
        ' {wildcards.thin}'
