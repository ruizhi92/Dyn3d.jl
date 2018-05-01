using FactCheck

facts("Testing add_body") do
    @fact add_body(nbody=4, ρ=0.01) --> body_system
