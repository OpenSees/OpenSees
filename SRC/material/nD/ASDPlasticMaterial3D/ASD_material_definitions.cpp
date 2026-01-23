

createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,YieldStress<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,YieldStress<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,YieldStress<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,YieldStress<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,YieldStress<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,YieldStress<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,YieldStress<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<TensorLinearHardeningFunction>, DP_cohesion<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,YieldStress<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<ArmstrongFrederickHardeningFunction>, DP_cohesion<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,YieldStress<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<TensorLinearHardeningFunction>, DP_cohesion<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,YieldStress<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<ArmstrongFrederickHardeningFunction>, DP_cohesion<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,YieldStress<ScalarLinearHardeningFunction>
            >, 
        MohrCoulomb_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,YieldStress<ScalarLinearHardeningFunction>
            >, 
        MohrCoulomb_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,DP_cohesion<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,DP_cohesion<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,DP_cohesion<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,DP_cohesion<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,DP_cohesion<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,DP_cohesion<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,DP_cohesion<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<TensorLinearHardeningFunction>, DP_cohesion<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,DP_cohesion<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<ArmstrongFrederickHardeningFunction>, DP_cohesion<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,DP_cohesion<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<TensorLinearHardeningFunction>, DP_cohesion<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,DP_cohesion<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<ArmstrongFrederickHardeningFunction>, DP_cohesion<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,DP_cohesion<ScalarLinearHardeningFunction>
            >, 
        MohrCoulomb_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,DP_cohesion<ScalarLinearHardeningFunction>
            >, 
        MohrCoulomb_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        MohrCoulomb_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        VonMises_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        MohrCoulomb_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        VonMises_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        MohrCoulomb_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        VonMises_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        MohrCoulomb_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        DruckerPrager_PF<
            BackStress<TensorLinearHardeningFunction>, DP_cohesion<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        MohrCoulomb_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        DruckerPrager_PF<
            BackStress<ArmstrongFrederickHardeningFunction>, DP_cohesion<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);



createASDPlasticMaterial3D<
        LinearIsotropic3D_EL, 
        MohrCoulomb_YF<
            BackStress<NullHardeningTensorFunction>
            >, 
        MohrCoulomb_PF<
            BackStress<NullHardeningTensorFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);

