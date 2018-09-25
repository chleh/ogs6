
# TODO [DUNE] FIX! APPEND_SOURCE_FILES(SOURCES)

add_library(BoundaryCondition
    BoundaryConditionCollection.cpp
    BoundaryConditionCollection.h
    BoundaryConditionConfig.h
    #BoundaryCondition.cpp
    BoundaryCondition.h
    #ConstraintDirichletBoundaryCondition.cpp
    #ConstraintDirichletBoundaryCondition.h
    #ConstraintDirichletBoundaryConditionLocalAssembler.h
    CreateBoundaryCondition.cpp
    CreateBoundaryCondition.h
    #DirichletBoundaryCondition.cpp
    #DirichletBoundaryCondition.h
    #GenericNaturalBoundaryCondition.h
    #GenericNaturalBoundaryCondition-impl.h
    #GenericNaturalBoundaryConditionLocalAssembler.h
    #GenericNonuniformNaturalBoundaryCondition.h
    #GenericNonuniformNaturalBoundaryCondition-impl.h
    #GenericNonuniformNaturalBoundaryConditionLocalAssembler.h
    #NeumannBoundaryCondition.cpp
    #NeumannBoundaryCondition.h
    #NeumannBoundaryConditionLocalAssembler.h
    #NonuniformDirichletBoundaryCondition.cpp
    NonuniformDirichletBoundaryConditionDUNE.cpp
    NonuniformDirichletBoundaryConditionDUNE.h
    #NonuniformDirichletBoundaryCondition.h
    #NonuniformNeumannBoundaryCondition.cpp
    #NonuniformNeumannBoundaryCondition.h
    #NonuniformNeumannBoundaryConditionLocalAssembler.h
    #NormalTractionBoundaryCondition.h
    #NormalTractionBoundaryCondition-impl.h
    #NormalTractionBoundaryConditionLocalAssembler.h
    #PhaseFieldIrreversibleDamageOracleBoundaryCondition.cpp
    #PhaseFieldIrreversibleDamageOracleBoundaryCondition.h
    #RobinBoundaryCondition.cpp
    #RobinBoundaryCondition.h
    #RobinBoundaryConditionLocalAssembler.h
    )

target_link_libraries(BoundaryCondition
    PUBLIC MeshLib NumLib)
