ni2, ni3, Ti2, Ti3, dlnnidr2, dlnnidr3, dlntidr2, dlntidr3 = TJLFEP.readEXPRO()

const exproConst::NamedTuple = (
    ni = (fill(NaN, 201), ni2, ni3, fill(NaN, 201)),
    Ti = (fill(NaN, 201), Ti2, Ti3, fill(NaN, 201)),
    dlnnidr = (fill(NaN, 201), dlnnidr2, dlnnidr3, fill(NaN, 201)),
    dlntidr = (fill(NaN, 201), dlntidr2, dlntidr3, fill(NaN, 201))
)