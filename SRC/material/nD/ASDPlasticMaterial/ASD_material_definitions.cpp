

createASDPlasticMaterial<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        LinearIsotropic3D_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        LinearIsotropic3D_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        LinearIsotropic3D_EL, 
        RoundedMohrCoulomb_YF<
            ScalarInternalVariable<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        LinearIsotropic3D_EL, 
        RoundedMohrCoulomb_YF<
            ScalarInternalVariable<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        LinearIsotropic3D_EL, 
        RoundedMohrCoulomb_YF<
            ScalarInternalVariable<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        LinearIsotropic3D_EL, 
        RoundedMohrCoulomb_YF<
            ScalarInternalVariable<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        LinearIsotropic3D_EL, 
        RoundedMohrCoulomb_YF<
            ScalarInternalVariable<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        LinearIsotropic3D_EL, 
        RoundedMohrCoulomb_YF<
            ScalarInternalVariable<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        DuncanChang_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        DuncanChang_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        DuncanChang_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        DuncanChang_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        DuncanChang_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        DuncanChang_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        DuncanChang_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        DuncanChang_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        DuncanChang_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        DuncanChang_EL, 
        VonMises_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        DuncanChang_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        DuncanChang_EL, 
        VonMises_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        DuncanChang_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        DuncanChang_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        DuncanChang_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        DuncanChang_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        DuncanChang_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        DuncanChang_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        DuncanChang_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        DuncanChang_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        DuncanChang_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        DuncanChang_EL, 
        DruckerPrager_YF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        DuncanChang_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        DuncanChang_EL, 
        DruckerPrager_YF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        DuncanChang_EL, 
        RoundedMohrCoulomb_YF<
            ScalarInternalVariable<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        DuncanChang_EL, 
        RoundedMohrCoulomb_YF<
            ScalarInternalVariable<ScalarLinearHardeningFunction>
            >, 
        VonMises_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        DuncanChang_EL, 
        RoundedMohrCoulomb_YF<
            ScalarInternalVariable<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<TensorLinearHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        DuncanChang_EL, 
        RoundedMohrCoulomb_YF<
            ScalarInternalVariable<ScalarLinearHardeningFunction>
            >, 
        DruckerPrager_PF<
            BackStress<ArmstrongFrederickHardeningFunction>,VonMisesRadius<ScalarLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        DuncanChang_EL, 
        RoundedMohrCoulomb_YF<
            ScalarInternalVariable<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<TensorLinearHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);




createASDPlasticMaterial<
        DuncanChang_EL, 
        RoundedMohrCoulomb_YF<
            ScalarInternalVariable<ScalarLinearHardeningFunction>
            >, 
        ConstantDilatancy_PF<
            BackStress<ArmstrongFrederickHardeningFunction>
            >
        > (instance_tag, yf_type, pf_type, el_type, iv_type, instance_pointers, available_models);


