ni1, ni2, ni3, ni4, Ti1, Ti2, Ti3, Ti4, dlnnidr1, dlnnidr2, dlnnidr3, dlnnidr4, dlntidr1, dlntidr2, dlntidr3, dlntidr4, cs, rmin_ex = TJLFEP.readEXPRO()

const exproConst::NamedTuple = (
    ni = (ni1, ni2, ni3, ni4),
    Ti = (Ti1, Ti2, Ti3, Ti4),
    dlnnidr = (dlnnidr1, dlnnidr2, dlnnidr3, dlnnidr4),
    dlntidr = (dlntidr1, dlntidr2, dlntidr3, dlntidr4),
    cs = cs,
    rmin_ex = rmin_ex
)