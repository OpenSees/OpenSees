#!/usr/bin/python
from itertools import product

EL = [
    "LinearIsotropic3D_EL",
    # "DuncanChang_EL",
]

YF = [
    "VonMises_YF",
    "DruckerPrager_YF",
    # "RoundedMohrCoulomb_YF",
]

PF = [
    "VonMises_PF",
    "DruckerPrager_PF",
    "ConstantDilatancy_PF"
]

IV_YF = {}

IV_YF["VonMises_YF"] = [
    "BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>",
    # "BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>",
]

IV_YF["DruckerPrager_YF"] = IV_YF["VonMises_YF"] 

IV_YF["RoundedMohrCoulomb_YF"] = [
    "ScalarInternalVariable<ScalarLinearHardeningFunction>"
]


#Options for PF variables depend on the model
IV_PF = {
    "VonMises_PF": [
    "BackStress<TensorLinearHardeningFunction>",
    "BackStress<ArmstrongFrederickHardeningFunction>",
    ],
    "DruckerPrager_PF":
    [
    "BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>",
    # "BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>",
    ]
}

IV_PF["ConstantDilatancy_PF"] = IV_PF["VonMises_PF"]





template = """

createASDPlasticMaterial<
        {EL}, 
        {YF}<
            {IV_YF}
            >, 
        {PF}<
            {IV_PF}
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);

"""


# with open("ASD_material_definitions.cpp","w") as fid:
#     for el, yf, pf, iv_yf, iv_pf in product(EL, YF, PF, IV_YF, IV_PF):
#         fid.write(template.format(EL=el, YF=yf, PF=pf, IV_YF=iv_yf, IV_PF=iv_pf))

with open("ASD_material_definitions.cpp","w") as fid:
    for el, yf, pf in product(EL, YF, PF):
        for iv_yf in IV_YF[yf]:
            for iv_pf in IV_PF[pf]:
                fid.write(template.format(EL=el, YF=yf, PF=pf, IV_YF=iv_yf, IV_PF=iv_pf))