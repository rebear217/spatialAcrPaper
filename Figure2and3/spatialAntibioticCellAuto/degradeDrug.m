function state = degradeDrug(state)
    state.drug = (1-state.drugDegradation)*state.drug;
    state.internalDrug = (1-state.drugDegradation)*state.internalDrug;
end