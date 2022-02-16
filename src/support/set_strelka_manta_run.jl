function set_strelka_manta_run(n_jo, me)

    ru = "--mode local --jobs $n_jo --memGb $me --quiet"

    return ru

end
