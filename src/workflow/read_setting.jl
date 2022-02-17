using JSON: parsefile

function read_setting(se)

    fe_va = parsefile(se)

    n_jo = fe_va["n_jo"]

    me = fe_va["me"]

    mo = fe_va["mo"]

    ta = fe_va["ta"]

    fr = fe_va["fr"]

    sd = fe_va["sd"]

    sa = fe_va["sa"]

    to = fe_va["to"]

    ou = fe_va["ou"]

    r1 = fe_va["r1"]

    r2 = fe_va["r2"]

    sor1 = fe_va["sor1"]

    sor2 = fe_va["sor2"]

    ge = fe_va["ge"]

    tr = fe_va["tr"]

    chs = fe_va["chs"]

    chn = fe_va["chn"]

    sn = fe_va["sn"]

    return n_jo, me, mo, ta, fr, sd, sa, to, ou, r1, r2, sor1, sor2, ge, tr, chs, chn, sn

end
