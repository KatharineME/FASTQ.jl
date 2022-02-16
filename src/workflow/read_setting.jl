using JSON: parsefile

function read_setting(se = joinpath(dirname(dirname(@__DIR__)), "input/setting_example.json"))

    fe_va = parsefile(se)

    n_jo = fe_va["n_jo"]

    me = fe_va["me"]

    mo = fe_va["mo"]

    ta = fe_va["ta"]

    sa = fe_va["sa"]

    to = fe_va["to"]

    ou = fe_va["ou"]

    ger1 = fe_va["ger1"]

    ger2 = fe_va["ger2"]

    sor1 = fe_va["sor1"]

    sor2 = fe_va["sor2"]

    ge = fe_va["ge"]

    tr = fe_va["tr"]

    chs = fe_va["chs"]

    chn = fe_va["chn"]

    sn = fe_va["sn"]

    return n_jo, me, mo, ta, sa, to, ou, ger1, ger2, sor1, sor2, ge, tr, chs, chn, sn

end
