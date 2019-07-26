"""Calculate the bearing load distribution on a lap joint

The lap joint is defined by a user-supplied *.csv file. The joint is assumed to
belayed out like the diagram below.


                          Lap 2
                 2-----4-----6-----8----> P2
                 |     |     |     |
                 B1    B2    B3    B4
                 |     |     |     |
        Pin <----1-----3-----5-----7----> P1
                          Lap 1

The bolts are shown as `Bi`, the nodes on the top lap are odd, and the
nodes on the lower lap are even. The input load is `Pin` and the output loads
are given by `P1` for the bottom lap and `P2` for the top lap.

Required modules:
    - math
    - csv
    - numpy
    - mechpy.bolt_flex

Inputs:


    The *.csv should be organized by fastener

# TODO: cleanup table samples
    fasteners
        ID    Thickness 1    Thickness 2    Modulus    Pitch    Diameter    Bolt Type    Nu
        1    t11    t12    e1    pitch1    d1    type1    nu_b_1
        2    t21    t22    e2    pitch2    d2    type2    nu_b_2
        3    t31    t32    e3    pitch3    d3    type3    nu_b_3
        4    t41    t42    e4    pitch4    d4    type4    nu_b_4

    joint info
        p_in        input load
        p_out1      lap 1 output load
        p_out2      lap 2 output load
        width       joint width
        e1          lap 1 modulus
        e2          lap 2 modulus
        mat1        lap 1 material type
        mat2        lap 2 material type
        nu1         lap 1 poisson's ratio
        nu2         lap2 poisson's ratio
        joint_type  single or double shear

# TODO: add complete list of formulations
    stiffness formulation
        - 'Huth'
        - etc.

Assumptions:
    - For stepped joints, the step occurs halfway between the fasteners


"""
