function state = replenishDrug(state)
    state.drug(state.drugLocationsX,state.drugLocationsY) = state.initialDrug;
end