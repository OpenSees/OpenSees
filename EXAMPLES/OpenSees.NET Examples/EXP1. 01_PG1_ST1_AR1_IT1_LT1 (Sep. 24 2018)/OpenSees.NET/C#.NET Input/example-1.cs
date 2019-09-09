using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace epfl_verifi1
{
    class Program
    {
        static void Main(string[] args)
        {
            var theDomain = new OpenSees.Components.DomainWrapper();

            #region paramaters

            var length = 1400.0;
            var width_Timoshenko_long = 1100.0;
            var width_Timoshenko_trans = 1400.0;
            var width_inside_long = 385.0;
            var width_inside_trans = 465.0;
            var width_Perimeter_long_right = 220.0;
            var width_Perimeter_long_left = 110.0;
            var width_Perimeter_trans_up = 250.0;
            var width_Perimeter_trans_down = 140.0;
            var tickness = 40.0;

            var Tension = 666.0;
            var IP = 8590.0;
            var OOP = 7177.0;
            var about_Tension = 1000000000000000.0;
            var about_IP = 1000000000000000.0; ;
            var about_OOP = 93000.0;

            var E_Timber_Long = 13200.0;
            var E_Timber_Trans = 2200.0;
            var poisson_ratio = 0.4;
            var G_Timber = 820.0;

            var Area_Timoshenko_long = width_Timoshenko_long * tickness;
            var Iz_Timoshenko_long = (tickness * width_Timoshenko_long * width_Timoshenko_long * width_Timoshenko_long) / 12.0;
            var Iy_Timoshenko_long = (width_Timoshenko_long * tickness * tickness * tickness) / (12.0);
            var Av_Timoshenko_long = (5 * width_Timoshenko_long * tickness) / (6);
            var Torsional_J_Timoshenko_long = (width_Timoshenko_long * tickness * (width_Timoshenko_long * width_Timoshenko_long + tickness * tickness)) / 12.0;

            var Area_Timoshenko_trans = width_Timoshenko_trans * tickness;
            var Iz_Timoshenko_trans = (tickness * width_Timoshenko_trans * width_Timoshenko_trans * width_Timoshenko_trans) / (12.0);
            var Iy_Timoshenko_trans = (width_Timoshenko_trans * tickness * tickness * tickness) / (12.0);
            var Av_Timoshenko_trans = (5 * width_Timoshenko_trans * tickness) / (6);
            var Torsional_J_Timoshenko_trans = (width_Timoshenko_trans * tickness * (width_Timoshenko_trans * width_Timoshenko_trans + tickness * tickness)) / 12.0;

            var Area_Perimeter_long_right = width_Perimeter_long_right * tickness;
            var Iz_Perimeter_long_right = 1.0 * (tickness * width_Perimeter_long_right * width_Perimeter_long_right * width_Perimeter_long_right) / (12.0);
            var Iy_Perimeter_long_right = 1.0 * (width_Perimeter_long_right * tickness * tickness * tickness) / (12.0);
            var Av_Perimeter_long_right = (5 * width_Perimeter_long_right * tickness) / (6);
            var Torsional_J_Perimeter_long_right = (width_Perimeter_long_right * tickness * (width_Perimeter_long_right * width_Perimeter_long_right + tickness * tickness)) / 12.0;

            var Area_Perimeter_long_left = width_Perimeter_long_left * tickness;
            var Iz_Perimeter_long_left = 1.0 * (tickness * width_Perimeter_long_left * width_Perimeter_long_left * width_Perimeter_long_left) / (12.0);
            var Iy_Perimeter_long_left = 1.0 * (width_Perimeter_long_left * tickness * tickness * tickness) / (12.0);
            var Av_Perimeter_long_left = (5 * width_Perimeter_long_left * tickness) / (6);
            var Torsional_J_Perimeter_long_left = (width_Perimeter_long_left * tickness * (width_Perimeter_long_left * width_Perimeter_long_left + tickness * tickness)) / 12.0;

            var Area_Perimeter_trans_up = width_Perimeter_trans_up * tickness;
            var Iz_Perimeter_trans_up = 1.0 * (tickness * width_Perimeter_trans_up * width_Perimeter_trans_up * width_Perimeter_trans_up) / (12.0);
            var Iy_Perimeter_trans_up = 1.0 * (width_Perimeter_trans_up * tickness * tickness * tickness) / (12.0);
            var Av_Perimeter_trans_up = (5 * width_Perimeter_trans_up * tickness) / (6);
            var Torsional_J_Perimeter_trans_up = (width_Perimeter_trans_up * tickness * (width_Perimeter_trans_up * width_Perimeter_trans_up + tickness * tickness)) / 12.0;

            var Area_Perimeter_trans_down = width_Perimeter_trans_down * tickness;
            var Iz_Perimeter_trans_down = 1.0 * (tickness * width_Perimeter_trans_down * width_Perimeter_trans_down * width_Perimeter_trans_down) / (12.0);
            var Iy_Perimeter_trans_down = 1.0 * (width_Perimeter_trans_down * tickness * tickness * tickness) / (12.0);
            var Av_Perimeter_trans_down = (5 * width_Perimeter_trans_down * tickness) / (6);
            var Torsional_J_Perimeter_trans_down = (width_Perimeter_trans_down * tickness * (width_Perimeter_trans_down * width_Perimeter_trans_down + tickness * tickness)) / 12.0;

            var Area_inside_long = width_inside_long * tickness;
            var Iz_inside_long = (tickness * width_inside_long * width_inside_long * width_inside_long) / (12.0);
            var Iy_inside_long = (width_inside_long * tickness * tickness * tickness) / (12.0);
            var Av_inside_long = (5 * width_inside_long * tickness) / (6);
            var Torsional_J_inside_long = (width_inside_long * tickness * (width_inside_long * width_inside_long + tickness * tickness)) / 12.0;

            var Area_inside_trans = width_inside_trans * tickness;
            var Iz_inside_trans = (tickness * width_inside_trans * width_inside_trans * width_inside_trans) / (12.0);
            var Iy_inside_trans = (width_inside_trans * tickness * tickness * tickness) / (12.0);
            var Av_inside_trans = (5 * width_inside_trans * tickness) / (6);
            var Torsional_J_inside_trans = (width_inside_trans * tickness * (width_inside_trans * width_inside_trans + tickness * tickness)) / 12.0;


            var rigidmattag = 1;
            var semi_rigidmattag = 101;

            var Freemattag = 2;

            var timbermattag_long = 30;
            var timbermattag_trans = 31;

            var Axial_perimeter_long_right = 40;
            var Axial_perimeter_long_left = 41;
            var Axial_perimeter_trans_up = 42;
            var Axial_perimeter_trans_down = 43;
            var Axial_inside_long = 44;
            var Axial_inside_trans = 45;

            var Shear_1_perimeter_long_right = 50;
            var Shear_1_perimeter_long_left = 51;
            var Shear_1_perimeter_trans_up = 52;
            var Shear_1_perimeter_trans_down = 53;
            var Shear_1_inside_long = 54;
            var Shear_1_inside_trans = 55;

            var Shear_2_perimeter_long = 60;
            var Shear_2_perimeter_trans = 61;
            var Shear_2_inside_long = 62;
            var Shear_2_inside_trans = 63;

            var Trans_Shear_Timoshenko_Eq = 70;
            var Long_Shear_Timoshenko_Eq = 71;

            var Stiffness_Tension = 10;
            var Stiffness_IP = 11;
            var Stiffness_OOP = 12;
            var Stiffness_about_Tension = 13;
            var Stiffness_about_IP = 14;
            var Stiffness_about_OOP = 15;
            #endregion

            #region materials

            var rigidmattagMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(rigidmattag, 1.0e+15, 0); //uniaxialMaterial Elastic $$rigidmattag[expr 1.0e+15]
            var semi_rigidmattagMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(semi_rigidmattag, 1.0e+4, 0);  //uniaxialMaterial Elastic $semi_rigidmattag[expr 1.0e+4]
            var FreemattagMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Freemattag, 1.0e-1, 0); // uniaxialMaterial Elastic $Freemattag[expr 1.0e-1]

            var timbermattag_longMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(timbermattag_long, E_Timber_Long, 0);//  uniaxialMaterial Elastic $timbermattag_long          $E_Timber_Long
            var timbermattag_transMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(timbermattag_trans, E_Timber_Trans, 0); //  uniaxialMaterial    Elastic $timbermattag_trans         $E_Timber_Trans
            var Axial_perimeter_long_rightMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Axial_perimeter_long_right, 100.0 * E_Timber_Long * Area_Perimeter_long_right / length, 0); // uniaxialMaterial    Elastic $Axial_perimeter_long_right[expr 100.0 *$E_Timber_Long *$Area_Perimeter_long_right /$length];         # for LT2
            var Axial_perimeter_long_leftMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Axial_perimeter_long_left, 100.0 * E_Timber_Long * Area_Perimeter_long_left / length, 0);// uniaxialMaterial Elastic $Axial_perimeter_long_left[expr 100.0 *$E_Timber_Long *$Area_Perimeter_long_left /$length];      # for LT2
            var Axial_perimeter_trans_upMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Axial_perimeter_trans_up, 100.0 * E_Timber_Trans * Area_Perimeter_trans_up / length, 0);// uniaxialMaterial Elastic $Axial_perimeter_trans_up[expr 100.0 *$E_Timber_Trans *$Area_Perimeter_trans_up /$length];      # for LT2
            var Axial_perimeter_trans_downMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Axial_perimeter_trans_down, 100.0 * E_Timber_Trans * Area_Perimeter_trans_down / length, 0);//uniaxialMaterial Elastic $Axial_perimeter_trans_down[expr 100.0 *$E_Timber_Trans *$Area_Perimeter_trans_down /$length];        # for LT2
            var Axial_inside_longMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Axial_inside_long, 1.0 * E_Timber_Long * Area_inside_long / length, 0); // uniaxialMaterial Elastic $Axial_inside_long[expr 1.*$E_Timber_Long *$Area_inside_long /$length];             # for LT2
            var Axial_inside_transMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Axial_inside_trans, 1.0 * E_Timber_Trans * Area_inside_trans / length, 0); //  uniaxialMaterial Elastic $Axial_inside_trans[expr 1.*$E_Timber_Trans *$Area_inside_trans /$length];           # for LT2


            var Shear_1_perimeter_long_rightMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Shear_1_perimeter_long_right, 100.0 * (5 * width_Perimeter_long_right * tickness) / (6), 0); // uniaxialMaterial Elastic $Shear_1_perimeter_long_right[expr 100.0 * (5 *$width_Perimeter_long_right *$tickness)/ (6)];
            var Shear_1_perimeter_long_leftMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Shear_1_perimeter_long_left, 100.0 * (5 * width_Perimeter_long_left * tickness) / (6), 0); // uniaxialMaterial Elastic $Shear_1_perimeter_long_left[expr 100.0 * (5 *$width_Perimeter_long_left *$tickness)/ (6)];
            var Shear_1_perimeter_trans_upMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Shear_1_perimeter_trans_up, 100.0 * (5 * width_Perimeter_trans_up * tickness) / (6), 0); // uniaxialMaterial Elastic $Shear_1_perimeter_trans_up[expr 100.0 * (5 *$width_Perimeter_trans_up *$tickness)/ (6)];

            var Shear_1_perimeter_trans_downMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Shear_1_perimeter_trans_down, 100.0 * (5 * width_Perimeter_trans_down * tickness) / (6), 0); // uniaxialMaterial Elastic $Shear_1_perimeter_trans_down[expr 100.0 * (5 *$width_Perimeter_trans_down *$tickness)/ (6)];
            var Shear_1_inside_longMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Shear_1_inside_long, 1.0 * (5 * width_inside_long * tickness) / (6), 0); // uniaxialMaterial Elastic $Shear_1_inside_long[expr 1.0 * (5 *$width_inside_long *$tickness)/ (6)];
            var Shear_1_inside_transMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Shear_1_inside_trans, 1.0 * (5 * width_inside_trans * tickness) / (6), 0); // uniaxialMaterial Elastic $Shear_1_inside_trans[expr 1.0 * (5 *$width_inside_trans *$tickness)/ (6)];

            var Shear_2_perimeter_longMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Shear_2_perimeter_long, 1.0e+15, 0);//uniaxialMaterial Elastic $Shear_2_perimeter_long     1.0e+15;
            var Shear_2_perimeter_transMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Shear_2_perimeter_trans, 1.0e+15, 0); // uniaxialMaterial Elastic $Shear_2_perimeter_trans    1.0e+15;
            var Shear_2_inside_longMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Shear_2_inside_long, 1.0e+15, 0); // uniaxialMaterial Elastic $Shear_2_inside_long        1.0e+15;
            var Shear_2_inside_transMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Shear_2_inside_trans, 1.0e+15, 0); // uniaxialMaterial Elastic $Shear_2_inside_trans       1.0e+15;


            var Trans_Shear_Timoshenko_EqMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Trans_Shear_Timoshenko_Eq, 34787.878, 0); // uniaxialMaterial Elastic $Trans_Shear_Timoshenko_Eq[expr 34787.878]
            var Long_Shear_Timoshenko_EqMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Long_Shear_Timoshenko_Eq, 16081.971, 0); // uniaxialMaterial Elastic $Long_Shear_Timoshenko_Eq[expr 16081.971]


            var Stiffness_TensionMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Stiffness_Tension, Tension, 0);//  uniaxialMaterial Elastic $Stiffness_Tension[expr   $Tension]
            var Stiffness_IPMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Stiffness_IP, IP, 0); // uniaxialMaterial Elastic $Stiffness_IP[expr   $IP]
            var Stiffness_OOPMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Stiffness_OOP, OOP, 0); // uniaxialMaterial Elastic $Stiffness_OOP[expr   $OOP]
            var Stiffness_about_TensionMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Stiffness_about_Tension, about_Tension, 0); // uniaxialMaterial Elastic $Stiffness_about_Tension[expr   $about_Tension]
            var Stiffness_about_IPMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Stiffness_about_IP, about_IP, 0); // uniaxialMaterial Elastic $Stiffness_about_IP[expr   $about_IP]
            var Stiffness_about_OOPMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Stiffness_about_OOP, about_OOP, 0); //uniaxialMaterial Elastic $Stiffness_about_OOP[expr   $about_OOP]

            #endregion


            #region nodes
            var node1 = new OpenSees.Components.NodeWrapper(1, 6, 1100, 0, 700);
            var node2 = new OpenSees.Components.NodeWrapper(2, 6, 1100, 0, 466.7);
            var node3 = new OpenSees.Components.NodeWrapper(3, 6, 1100, 0, 765);
            var node4 = new OpenSees.Components.NodeWrapper(4, 6, 0, 0, 700);
            var node5 = new OpenSees.Components.NodeWrapper(5, 6, 0, 0, 933.4);
            var node6 = new OpenSees.Components.NodeWrapper(6, 6, 0, 0, 655);
            var node7 = new OpenSees.Components.NodeWrapper(7, 6, 1100, 0, 1155);
            var node8 = new OpenSees.Components.NodeWrapper(8, 6, 1100, 0, 1400);
            var node9 = new OpenSees.Components.NodeWrapper(9, 6, 1100, 0, 933.4);
            var node10 = new OpenSees.Components.NodeWrapper(10, 6, 1100, 0, 355);
            var node11 = new OpenSees.Components.NodeWrapper(11, 6, 1100, 0, 0);
            var node12 = new OpenSees.Components.NodeWrapper(12, 6, 0, 0, 1045);
            var node13 = new OpenSees.Components.NodeWrapper(13, 6, 0, 0, 1400);
            var node14 = new OpenSees.Components.NodeWrapper(14, 6, 0, 0, 466.7);
            var node15 = new OpenSees.Components.NodeWrapper(15, 6, 0, 0, 245);
            var node16 = new OpenSees.Components.NodeWrapper(16, 6, 0, 0, 0);
            var node17 = new OpenSees.Components.NodeWrapper(17, 6, 1140, 0, 355);
            var node18 = new OpenSees.Components.NodeWrapper(18, 6, 1140, 0, 765);
            var node19 = new OpenSees.Components.NodeWrapper(19, 6, 1140, 0, 1155);
            var node20 = new OpenSees.Components.NodeWrapper(20, 6, 255.207, 0, 1440);
            var node21 = new OpenSees.Components.NodeWrapper(21, 6, 255.207, 0, 1400);
            var node22 = new OpenSees.Components.NodeWrapper(22, 6, 672.394, 0, 1440);
            var node23 = new OpenSees.Components.NodeWrapper(23, 6, 672.394, 0, 1400);
            var node24 = new OpenSees.Components.NodeWrapper(24, 6, -40, 0, 1045);
            var node25 = new OpenSees.Components.NodeWrapper(25, 6, -40, 0, 655);
            var node26 = new OpenSees.Components.NodeWrapper(26, 6, -40, 0, 245);
            var node27 = new OpenSees.Components.NodeWrapper(27, 6, 426.25, 0, -40);
            var node28 = new OpenSees.Components.NodeWrapper(28, 6, 426.25, 0, 0);
            var node29 = new OpenSees.Components.NodeWrapper(29, 6, 838.75, 0, -40);
            var node30 = new OpenSees.Components.NodeWrapper(30, 6, 838.75, 0, 0);
            var node31 = new OpenSees.Components.NodeWrapper(31, 6, 733.34, 0, 1400);
            var node32 = new OpenSees.Components.NodeWrapper(32, 6, 366.67, 0, 1400);
            var node33 = new OpenSees.Components.NodeWrapper(33, 6, 550, 0, 1400);
            var node34 = new OpenSees.Components.NodeWrapper(34, 6, 733.34, 0, 0);
            var node35 = new OpenSees.Components.NodeWrapper(35, 6, 550, 0, 0);
            var node36 = new OpenSees.Components.NodeWrapper(36, 6, 366.67, 0, 0);
            var node1601 = new OpenSees.Components.NodeWrapper(1601, 6, 0, 0, 0);
            var node1602 = new OpenSees.Components.NodeWrapper(1602, 6, 0, 0, 0);
            var node1301 = new OpenSees.Components.NodeWrapper(1301, 6, 0, 0, 1400);
            var node1302 = new OpenSees.Components.NodeWrapper(1302, 6, 0, 0, 1400);
            var node1101 = new OpenSees.Components.NodeWrapper(1101, 6, 1100, 0, 0);
            var node1102 = new OpenSees.Components.NodeWrapper(1102, 6, 1100, 0, 0);
            var node802 = new OpenSees.Components.NodeWrapper(802, 6, 1100, 0, 1400);
            var node801 = new OpenSees.Components.NodeWrapper(801, 6, 1100, 0, 1400);
            var node3201 = new OpenSees.Components.NodeWrapper(3201, 6, 366.67, 0, 1400);
            var node3101 = new OpenSees.Components.NodeWrapper(3101, 6, 733.34, 0, 1400);
            var node3601 = new OpenSees.Components.NodeWrapper(3601, 6, 366.67, 0, 0);
            var node3401 = new OpenSees.Components.NodeWrapper(3401, 6, 733.34, 0, 0);
            var node1402 = new OpenSees.Components.NodeWrapper(1402, 6, 0, 0, 466.7);
            var node502 = new OpenSees.Components.NodeWrapper(502, 6, 0, 0, 933.4);
            var node902 = new OpenSees.Components.NodeWrapper(902, 6, 1100, 0, 933.4);
            var node202 = new OpenSees.Components.NodeWrapper(202, 6, 1100, 0, 466.7);

            var nodes = new OpenSees.Components.NodeWrapper[] { node1, node2, node3, node4, node5, node6, node7, node8, node9, node10,
            node11, node12, node13, node14, node15, node16, node17, node18, node19, node20,
            node21, node22, node23, node24, node25, node26, node27, node28, node29, node30,
            node31, node32, node33, node34, node35, node36 ,
            node1601, node1602, node1301, node1302, node1101, node1102, node802, node801, node3201, node3101,
            node3601, node3401, node1402, node502, node902, node202};
            theDomain.AddNode(nodes);

            foreach (var tag in new[] { 27, 29 /*, 24,25,26*/ })
                for (var dof = 0; dof < 6; dof++)
                    theDomain.AddSP_Constraint(new OpenSees.Components.Constraints.SP_ConstraintWrapper(tag, dof, 0.0, true));


            #endregion

            #region elements
            var xzplaneVec = new OpenSees.VectorWrapper(new double[] { 0, 1, 0 });
            var lineartag_TimoshenkoCrdT = new OpenSees.Elements.CrdTransfs.LinearCrdTransf3dWrapper(xzplaneVec);
            var lineartag_VerticalCrdT = new OpenSees.Elements.CrdTransfs.LinearCrdTransf3dWrapper(xzplaneVec);
            var lineartag_HorizontalCrdT = new OpenSees.Elements.CrdTransfs.LinearCrdTransf3dWrapper(xzplaneVec);

            {
                var dirIds1 = new OpenSees.IDWrapper(new int[] { 1 }); // 2-- for wrapper
                var ele_link_tnl1 = new OpenSees.Elements.TwoNodeLinkWrapper(1, 3, 35, 33, dirIds1, new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Long_Shear_Timoshenko_EqMat }, new OpenSees.VectorWrapper(new double[] { 1, 0, 0 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                var ele_link_tnl58 = new OpenSees.Elements.TwoNodeLinkWrapper(58, 3, 1, 4, dirIds1, new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Trans_Shear_Timoshenko_EqMat }, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(ele_link_tnl1);
                theDomain.AddElement(ele_link_tnl58);
            }
            

            //puts "The Timoshenko beams are defined"
            //# Defining Link Elements Between Zero-Length Nodes
            var dirIds2 = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });


            {
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Axial_perimeter_long_leftMat, Shear_1_perimeter_long_leftMat, rigidmattagMat, FreemattagMat, FreemattagMat, FreemattagMat };
                var ele_link_tnl2 = new OpenSees.Elements.TwoNodeLinkWrapper(2, 3, 16, 1601, dirIds2, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 1, 0, 0 }), new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(ele_link_tnl2);
            }

            {
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Axial_perimeter_trans_downMat, Shear_1_perimeter_trans_downMat, rigidmattagMat, FreemattagMat, FreemattagMat, FreemattagMat };
                var ele_link_tnl3 = new OpenSees.Elements.TwoNodeLinkWrapper(3, 3, 1602, 16, dirIds2, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(new double[] { -1, 0, 0 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(ele_link_tnl3);
            }

            {
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Axial_perimeter_long_leftMat, Shear_1_perimeter_long_leftMat, rigidmattagMat, FreemattagMat, FreemattagMat, FreemattagMat };
                var ele_link_tnl4 = new OpenSees.Elements.TwoNodeLinkWrapper(4, 3, 1301, 13, dirIds2, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 1, 0, 0 }), new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(ele_link_tnl4);
            }

            {
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Axial_perimeter_trans_upMat, Shear_1_perimeter_trans_upMat, rigidmattagMat, FreemattagMat, FreemattagMat, FreemattagMat };
                var ele_link_tnl5 = new OpenSees.Elements.TwoNodeLinkWrapper(5, 3, 1302, 13, dirIds2, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(new double[] { -1, 0, 0 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(ele_link_tnl5);
            }

            {
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Axial_perimeter_long_rightMat, Shear_1_perimeter_long_rightMat, rigidmattagMat, FreemattagMat, FreemattagMat, FreemattagMat };
                var ele_link_tnl6 = new OpenSees.Elements.TwoNodeLinkWrapper(6, 3, 11, 1101, dirIds2, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 1, 0, 0 }), new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(ele_link_tnl6);
            }

            {
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Axial_perimeter_trans_downMat, Shear_1_perimeter_trans_downMat, rigidmattagMat, FreemattagMat, FreemattagMat, FreemattagMat };
                var ele_link_tnl7 = new OpenSees.Elements.TwoNodeLinkWrapper(7, 3, 11, 1102, dirIds2, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(new double[] { -1, 0, 0 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(ele_link_tnl7);
            }

            {
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Axial_perimeter_long_rightMat, Shear_1_perimeter_long_rightMat, rigidmattagMat, FreemattagMat, FreemattagMat, FreemattagMat };
                var ele_link_tnl8 = new OpenSees.Elements.TwoNodeLinkWrapper(8, 3, 801, 8, dirIds2, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 1, 0,0 }), new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(ele_link_tnl8);
            }

            {
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Axial_perimeter_trans_upMat, Shear_1_perimeter_trans_upMat, rigidmattagMat, FreemattagMat, FreemattagMat, FreemattagMat };
                var ele_link_tnl9 = new OpenSees.Elements.TwoNodeLinkWrapper(9, 3, 8, 802, dirIds2, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] {0,0,1 }), new OpenSees.VectorWrapper(new double[] { -1,0,0 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(ele_link_tnl9);
            }

            {
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Axial_inside_transMat, Shear_1_inside_transMat, rigidmattagMat, FreemattagMat, FreemattagMat, FreemattagMat };
                var ele_link_tnl10 = new OpenSees.Elements.TwoNodeLinkWrapper(10, 3, 2, 202, dirIds2, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(new double[] { -1, 0, 0 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                var ele_link_tnl11 = new OpenSees.Elements.TwoNodeLinkWrapper(11, 3, 9, 902, dirIds2, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(new double[] { -1, 0, 0 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                var ele_link_tnl12 = new OpenSees.Elements.TwoNodeLinkWrapper(12, 3, 502, 5, dirIds2, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(new double[] { -1, 0, 0 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                var ele_link_tnl13 = new OpenSees.Elements.TwoNodeLinkWrapper(13, 3, 1402, 14, dirIds2, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(new double[] { -1, 0, 0 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(ele_link_tnl10);
                theDomain.AddElement(ele_link_tnl11);
                theDomain.AddElement(ele_link_tnl12);
                theDomain.AddElement(ele_link_tnl13);
            }

            {
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Axial_inside_longMat, Shear_1_inside_longMat, rigidmattagMat, FreemattagMat, FreemattagMat, FreemattagMat };
                var ele_link_tnl14 = new OpenSees.Elements.TwoNodeLinkWrapper(14, 3, 36, 3601, dirIds2, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 1, 0, 0 }), new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                var ele_link_tnl15 = new OpenSees.Elements.TwoNodeLinkWrapper(15, 3, 34, 3401, dirIds2, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 1, 0, 0 }), new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                var ele_link_tnl16 = new OpenSees.Elements.TwoNodeLinkWrapper(16, 3, 3201, 32, dirIds2, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 1, 0, 0 }), new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                var ele_link_tnl17 = new OpenSees.Elements.TwoNodeLinkWrapper(17, 3, 3101, 31, dirIds2, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 1, 0, 0 }), new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(ele_link_tnl14);
                theDomain.AddElement(ele_link_tnl15);
                theDomain.AddElement(ele_link_tnl16);
                theDomain.AddElement(ele_link_tnl17);
            }

            //# Defining Vertical Perimeter Elements LEFT
            {
                var ele_timbeam18 = new OpenSees.Elements.ElasticTimoshenkoBeam3dWrapper(18, 1601, 15, E_Timber_Long, G_Timber, Area_Perimeter_long_left, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left, Av_Perimeter_long_left, Av_Perimeter_long_left, lineartag_VerticalCrdT, 0, 0);
                var ele_timbeam19 = new OpenSees.Elements.ElasticTimoshenkoBeam3dWrapper(19, 15, 14, E_Timber_Long, G_Timber, Area_Perimeter_long_left, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left, Av_Perimeter_long_left, Av_Perimeter_long_left, lineartag_VerticalCrdT, 0, 0);
                var ele_timbeam20 = new OpenSees.Elements.ElasticTimoshenkoBeam3dWrapper(20, 14, 6, E_Timber_Long, G_Timber, Area_Perimeter_long_left, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left, Av_Perimeter_long_left, Av_Perimeter_long_left, lineartag_VerticalCrdT, 0, 0);
                var ele_timbeam21 = new OpenSees.Elements.ElasticTimoshenkoBeam3dWrapper(21, 6, 4, E_Timber_Long, G_Timber, Area_Perimeter_long_left, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left, Av_Perimeter_long_left, Av_Perimeter_long_left, lineartag_VerticalCrdT, 0, 0);
                var ele_timbeam22 = new OpenSees.Elements.ElasticTimoshenkoBeam3dWrapper(22, 4, 5, E_Timber_Long, G_Timber, Area_Perimeter_long_left, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left, Av_Perimeter_long_left, Av_Perimeter_long_left, lineartag_VerticalCrdT, 0, 0);
                var ele_timbeam23 = new OpenSees.Elements.ElasticTimoshenkoBeam3dWrapper(23, 5, 12, E_Timber_Long, G_Timber, Area_Perimeter_long_left, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left, Av_Perimeter_long_left, Av_Perimeter_long_left, lineartag_VerticalCrdT, 0, 0);
                var ele_timbeam24 = new OpenSees.Elements.ElasticTimoshenkoBeam3dWrapper(24, 12, 1301, E_Timber_Long, G_Timber, Area_Perimeter_long_left, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left, Av_Perimeter_long_left, Av_Perimeter_long_left, lineartag_VerticalCrdT, 0, 0);
                theDomain.AddElement(ele_timbeam18);
                theDomain.AddElement(ele_timbeam19);
                theDomain.AddElement(ele_timbeam20);
                theDomain.AddElement(ele_timbeam21);
                theDomain.AddElement(ele_timbeam22);
                theDomain.AddElement(ele_timbeam23);
                theDomain.AddElement(ele_timbeam24);
            }

            //# Defining Vertical Perimeter Elements RIGHT
            {
                var ele_timbeam25 = new OpenSees.Elements.ElasticTimoshenkoBeam3dWrapper(25, 1101, 10, E_Timber_Long, G_Timber, Area_Perimeter_long_right, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right, Av_Perimeter_long_right, Av_Perimeter_long_right, lineartag_VerticalCrdT, 0, 0);
                var ele_timbeam26 = new OpenSees.Elements.ElasticTimoshenkoBeam3dWrapper(26, 10, 2, E_Timber_Long, G_Timber, Area_Perimeter_long_right, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right, Av_Perimeter_long_right, Av_Perimeter_long_right, lineartag_VerticalCrdT, 0, 0);
                var ele_timbeam27 = new OpenSees.Elements.ElasticTimoshenkoBeam3dWrapper(27, 2, 1, E_Timber_Long, G_Timber, Area_Perimeter_long_right, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right, Av_Perimeter_long_right, Av_Perimeter_long_right, lineartag_VerticalCrdT, 0, 0);
                var ele_timbeam28 = new OpenSees.Elements.ElasticTimoshenkoBeam3dWrapper(28, 1, 3, E_Timber_Long, G_Timber, Area_Perimeter_long_right, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right, Av_Perimeter_long_right, Av_Perimeter_long_right, lineartag_VerticalCrdT, 0, 0);
                var ele_timbeam29 = new OpenSees.Elements.ElasticTimoshenkoBeam3dWrapper(29, 3, 9, E_Timber_Long, G_Timber, Area_Perimeter_long_right, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right, Av_Perimeter_long_right, Av_Perimeter_long_right, lineartag_VerticalCrdT, 0, 0);
                var ele_timbeam30 = new OpenSees.Elements.ElasticTimoshenkoBeam3dWrapper(30, 9, 7, E_Timber_Long, G_Timber, Area_Perimeter_long_right, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right, Av_Perimeter_long_right, Av_Perimeter_long_right, lineartag_VerticalCrdT, 0, 0);
                var ele_timbeam31 = new OpenSees.Elements.ElasticTimoshenkoBeam3dWrapper(31, 7, 801, E_Timber_Long, G_Timber, Area_Perimeter_long_right, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right, Av_Perimeter_long_right, Av_Perimeter_long_right, lineartag_VerticalCrdT, 0, 0);
                theDomain.AddElement(ele_timbeam25);
                theDomain.AddElement(ele_timbeam26);
                theDomain.AddElement(ele_timbeam27);
                theDomain.AddElement(ele_timbeam28);
                theDomain.AddElement(ele_timbeam29);
                theDomain.AddElement(ele_timbeam30);
                theDomain.AddElement(ele_timbeam31);
            }

            //# Defining Vertical Perimeter Elements DOWN
            {
                var ele_timbeam32 = new OpenSees.Elements.ElasticTimoshenkoBeam3dWrapper(32, 1102, 30, E_Timber_Trans, G_Timber, Area_Perimeter_trans_down, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down, Av_Perimeter_trans_down, Av_Perimeter_trans_down, lineartag_HorizontalCrdT, 0, 0);
                var ele_timbeam33 = new OpenSees.Elements.ElasticTimoshenkoBeam3dWrapper(33, 30, 34, E_Timber_Trans, G_Timber, Area_Perimeter_trans_down, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down, Av_Perimeter_trans_down, Av_Perimeter_trans_down, lineartag_HorizontalCrdT, 0, 0);
                var ele_timbeam34 = new OpenSees.Elements.ElasticTimoshenkoBeam3dWrapper(34, 34, 35, E_Timber_Trans, G_Timber, Area_Perimeter_trans_down, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down, Av_Perimeter_trans_down, Av_Perimeter_trans_down, lineartag_HorizontalCrdT, 0, 0);
                var ele_timbeam35 = new OpenSees.Elements.ElasticTimoshenkoBeam3dWrapper(35, 35, 28, E_Timber_Trans, G_Timber, Area_Perimeter_trans_down, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down, Av_Perimeter_trans_down, Av_Perimeter_trans_down, lineartag_HorizontalCrdT, 0, 0);
                var ele_timbeam36 = new OpenSees.Elements.ElasticTimoshenkoBeam3dWrapper(36, 28, 36, E_Timber_Trans, G_Timber, Area_Perimeter_trans_down, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down, Av_Perimeter_trans_down, Av_Perimeter_trans_down, lineartag_HorizontalCrdT, 0, 0);
                var ele_timbeam37 = new OpenSees.Elements.ElasticTimoshenkoBeam3dWrapper(37, 36, 1602, E_Timber_Trans, G_Timber, Area_Perimeter_trans_down, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down, Av_Perimeter_trans_down, Av_Perimeter_trans_down, lineartag_HorizontalCrdT, 0, 0);
                
                theDomain.AddElement(ele_timbeam32);
                theDomain.AddElement(ele_timbeam33);
                theDomain.AddElement(ele_timbeam34);
                theDomain.AddElement(ele_timbeam35);
                theDomain.AddElement(ele_timbeam36);
                theDomain.AddElement(ele_timbeam37);
                
            }

            //# Defining Vertical Perimeter Elements UP
            {
                var ele_timbeam38 = new OpenSees.Elements.ElasticTimoshenkoBeam3dWrapper(38, 802, 31, E_Timber_Trans, G_Timber, Area_Perimeter_trans_up, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up, Av_Perimeter_trans_up, Av_Perimeter_trans_up, lineartag_HorizontalCrdT, 0, 0);
                var ele_timbeam39 = new OpenSees.Elements.ElasticTimoshenkoBeam3dWrapper(39, 31, 23, E_Timber_Trans, G_Timber, Area_Perimeter_trans_up, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up, Av_Perimeter_trans_up, Av_Perimeter_trans_up, lineartag_HorizontalCrdT, 0, 0);
                var ele_timbeam40 = new OpenSees.Elements.ElasticTimoshenkoBeam3dWrapper(40, 23, 33, E_Timber_Trans, G_Timber, Area_Perimeter_trans_up, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up, Av_Perimeter_trans_up, Av_Perimeter_trans_up, lineartag_HorizontalCrdT, 0, 0);
                var ele_timbeam41 = new OpenSees.Elements.ElasticTimoshenkoBeam3dWrapper(41, 33, 32, E_Timber_Trans, G_Timber, Area_Perimeter_trans_up, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up, Av_Perimeter_trans_up, Av_Perimeter_trans_up, lineartag_HorizontalCrdT, 0, 0);
                var ele_timbeam42 = new OpenSees.Elements.ElasticTimoshenkoBeam3dWrapper(42, 32, 21, E_Timber_Trans, G_Timber, Area_Perimeter_trans_up, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up, Av_Perimeter_trans_up, Av_Perimeter_trans_up, lineartag_HorizontalCrdT, 0, 0);
                var ele_timbeam43 = new OpenSees.Elements.ElasticTimoshenkoBeam3dWrapper(43, 21, 1302, E_Timber_Trans, G_Timber, Area_Perimeter_trans_up, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up, Av_Perimeter_trans_up, Av_Perimeter_trans_up, lineartag_HorizontalCrdT, 0, 0);

                theDomain.AddElement(ele_timbeam38);
                theDomain.AddElement(ele_timbeam39);
                theDomain.AddElement(ele_timbeam40);
                theDomain.AddElement(ele_timbeam41);
                theDomain.AddElement(ele_timbeam42);
                theDomain.AddElement(ele_timbeam43);
            }

            //# Defining Vertical Inside Elements
            {
                var ele_timbeam44 = new OpenSees.Elements.ElasticTimoshenkoBeam3dWrapper(44, 3601, 3201, E_Timber_Long, G_Timber, Area_inside_long, Torsional_J_inside_long, Iy_inside_long, Iz_inside_long, Av_inside_long, Av_inside_long, lineartag_VerticalCrdT, 0, 0);
                var ele_timbeam45 = new OpenSees.Elements.ElasticTimoshenkoBeam3dWrapper(45, 3401, 3101, E_Timber_Long, G_Timber, Area_inside_long, Torsional_J_inside_long, Iy_inside_long, Iz_inside_long, Av_inside_long, Av_inside_long, lineartag_VerticalCrdT, 0, 0);

                theDomain.AddElement(ele_timbeam44);
                theDomain.AddElement(ele_timbeam45);
            }

            //# Defining Vertical Inside Elements
            {
                var ele_timbeam46 = new OpenSees.Elements.ElasticTimoshenkoBeam3dWrapper(46, 202, 1402, E_Timber_Trans, G_Timber, Area_inside_trans, Torsional_J_inside_trans, Iy_inside_trans, Iz_inside_trans, Av_inside_trans, Av_inside_trans, lineartag_HorizontalCrdT, 0, 0);
                var ele_timbeam47 = new OpenSees.Elements.ElasticTimoshenkoBeam3dWrapper(47, 902, 502, E_Timber_Trans, G_Timber, Area_inside_trans, Torsional_J_inside_trans, Iy_inside_trans, Iz_inside_trans, Av_inside_trans, Av_inside_trans, lineartag_HorizontalCrdT, 0, 0);

                theDomain.AddElement(ele_timbeam46);
                theDomain.AddElement(ele_timbeam47);
            }

            //# Defining Joint Elements
            {
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var ele_link_tnl48 = new OpenSees.Elements.TwoNodeLinkWrapper(48, 3, 27, 28, dirIds2, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 1, 0, 0 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                var ele_link_tnl49 = new OpenSees.Elements.TwoNodeLinkWrapper(49, 3, 29, 30, dirIds2, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 1, 0, 0 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                var ele_link_tnl50 = new OpenSees.Elements.TwoNodeLinkWrapper(50, 3, 23, 22, dirIds2, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 1, 0, 0 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                var ele_link_tnl51 = new OpenSees.Elements.TwoNodeLinkWrapper(51, 3, 21, 20, dirIds2, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 1, 0, 0 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(ele_link_tnl48);
                theDomain.AddElement(ele_link_tnl49);
                theDomain.AddElement(ele_link_tnl50);
                theDomain.AddElement(ele_link_tnl51);
            }

            {
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var ele_link_tnl52 = new OpenSees.Elements.TwoNodeLinkWrapper(52, 3, 19, 7, dirIds2, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                var ele_link_tnl53 = new OpenSees.Elements.TwoNodeLinkWrapper(53, 3, 18, 3, dirIds2, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                var ele_link_tnl54 = new OpenSees.Elements.TwoNodeLinkWrapper(54, 3, 17, 10, dirIds2, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);

                var ele_link_tnl55 = new OpenSees.Elements.TwoNodeLinkWrapper(55, 3, 12, 24, dirIds2, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                var ele_link_tnl56 = new OpenSees.Elements.TwoNodeLinkWrapper(56, 3, 6, 25, dirIds2, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                var ele_link_tnl57 = new OpenSees.Elements.TwoNodeLinkWrapper(57, 3, 15, 26, dirIds2, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);

                theDomain.AddElement(ele_link_tnl52);
                theDomain.AddElement(ele_link_tnl53);
                theDomain.AddElement(ele_link_tnl54);

                theDomain.AddElement(ele_link_tnl55);
                theDomain.AddElement(ele_link_tnl56);
                theDomain.AddElement(ele_link_tnl57);

            }

            #endregion

            #region recorders
            var savepath = System.Environment.CurrentDirectory + @"\opsnet_results\";
            if (!System.IO.Directory.Exists(savepath))
                System.IO.Directory.CreateDirectory(savepath);


            {
                var nodeTags = new int[] { 21, 23, 33, 13, 8, 7, 3, 10, 19, 18, 17, 20, 22, 28 };
                foreach(var tag in nodeTags)
                {
                    var opsstream = new OpenSees.Handlers.DataFileStreamWrapper(savepath + $@"\OPSNET_Disp_{tag}.txt");
                    var recorder = new OpenSees.Recorders.NodeRecorderWrapper(dirIds2, new OpenSees.IDWrapper(new int[] { tag }), 0, "disp", theDomain, opsstream);
                    theDomain.AddRecorder(recorder);
                }
            }

            {
                var nodeTags = new int[] { 27,29 };
                foreach(var tag in nodeTags)
                {
                    var opsstream = new OpenSees.Handlers.DataFileStreamWrapper(savepath + $@"\OPSNET_Nodal_Reactions_{tag}.txt");
                    var recorder = new OpenSees.Recorders.NodeRecorderWrapper(dirIds2, new OpenSees.IDWrapper(new int[] { tag }), 0, "reaction", theDomain, opsstream);
                    theDomain.AddRecorder(recorder);
                }
                
            }

            {
                var eleTags = new int[] { 1, 2 };

                foreach (var tag in eleTags)
                {
                    var opsstream = new OpenSees.Handlers.DataFileStreamWrapper(savepath + $@"\OPSNET_Element_Forces_{tag}.txt");
                    var recorder = new OpenSees.Recorders.ElementRecorderWrapper(new OpenSees.IDWrapper(new int[] { tag }), new string[] { "force" }, true, theDomain, opsstream, 0, new OpenSees.IDWrapper(0));
                    theDomain.AddRecorder(recorder);
                }
                
            }
            #endregion

            #region loading 
            var theSeries = new OpenSees.Components.Timeseries.LinearSeriesWrapper();
            var theLoadPattern = new OpenSees.Components.LoadPatterns.LoadPatternWrapper(1);
            theLoadPattern.SetTimeSeries(theSeries);
            theDomain.AddLoadPattern(theLoadPattern);

            var theLoadValues = new OpenSees.VectorWrapper(6);
            theLoadValues.Zero();
            theLoadValues[0] = 9871.8;

            theDomain.AddNodalLoad(new OpenSees.Components.Loads.NodalLoadWrapper(1,20, theLoadValues, false), 1);
            theDomain.AddNodalLoad(new OpenSees.Components.Loads.NodalLoadWrapper(2,22, theLoadValues, false), 1);

            var theModel = new OpenSees.AnalysisModelWrapper();
            var theSolnAlgo = new OpenSees.Algorithms.NewtonRaphsonWrapper();
            var theIntegrator = new OpenSees.Integrators.Static.LoadControlWrapper(0.0005, 1, 0.0005, 0.0005);
            var theHandler = new OpenSees.Handlers.PlainHandlerWrapper();
            var theRCM = new OpenSees.GraphNumberers.RCMWrapper(false);
            var theNumberer = new OpenSees.Numberers.DOF_NumbererWrapper(theRCM);
            var theSolver = new OpenSees.Systems.Linears.BandGenLinLapackSolverWrapper();
            var theSOE = new OpenSees.Systems.Linears.BandGenLinSOEWrapper(theSolver);
            var test = new OpenSees.ConvergenceTests.CTestNormDispIncrWrapper(1e-12, 20, 2, 2, 1.0e10);
            var theAnalysis = new OpenSees.Analysis.StaticAnalysisWrapper(
                theDomain,
                theHandler,
                theNumberer,
                theModel,
                theSolnAlgo,
                theSOE,
                theIntegrator,
                test);

            theAnalysis.Analyze(2000);
            

            Console.ReadKey();
            #endregion
        }
    }
}
