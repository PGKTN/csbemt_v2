using System;
using System.Collections.Generic;
using System.Dynamic;
using System.Linq;
using System.Text;
using System.Text.Encodings.Web;
using System.Threading.Tasks;

namespace csbemt_v2
{
    internal class Program
    {
        static void Main()
        {
            ReadWriteINIfile ini = new ReadWriteINIfile(@"C:\csbemt_v2\rotor.ini");

            double rho = double.Parse(ini.ReadINI("fluid", "rho"));
            double mu = double.Parse(ini.ReadINI("fluid", "mu"));

            double rpm = double.Parse(ini.ReadINI("case", "rpm"));
            double v_inf = double.Parse(ini.ReadINI("case", "v_inf"));

            string[] sections = ini.ReadINI("rotor", "section").Split(' ');

            string[] str_radius = ini.ReadINI("rotor", "radius").Split(' ');
            string[] str_chord = ini.ReadINI("rotor", "chord").Split(' ');
            string[] str_theta = ini.ReadINI("rotor", "pitch").Split(' ');

            int Nb = int.Parse(ini.ReadINI("rotor", "nblades"));

            List<double> radius = new List<double>();
            List<double> chord = new List<double>();
            List<double> theta = new List<double>();

            double omega = 2 * Math.PI * rpm / 60;
            Console.WriteLine("omega : " + omega);

            Console.WriteLine("rho : " + rho);
            Console.WriteLine("Nb : " + Nb);

            Calculate Calc = new Calculate();

            for (int i = 0; i < sections.Length; i++)
            {
                radius.Add(double.Parse(str_radius[i]));
                chord.Add(double.Parse(str_chord[i]));
                theta.Add(double.Parse(str_theta[i]));
            }

            for (int i = 0; i < sections.Length; i++)
                Console.WriteLine("radius[" + i + "] : " + radius[i]);
            Console.WriteLine();

            for (int i = 0; i < sections.Length; i++)
                Console.WriteLine("chord[" + i + "] : " + chord[i]);
            Console.WriteLine();

            for (int i = 0; i < sections.Length; i++)
                Console.WriteLine("twist[" + i + "] : " + theta[i]);

            List<double> reynolds = new List<double>();
            for (int i = 0; i < sections.Length; i++)
            {
                reynolds.Add(Calc.Get_Reynolds(rho, omega, radius[i], chord[i], mu));
                Console.WriteLine("Reynolds[" + i + "] : " + reynolds[i]);
            }

            Console.WriteLine();
            double sigma = 0.0;
            sigma = (Nb * chord[sections.Length - 1] * (radius[sections.Length - 1] - radius[0])) / (Math.PI * Math.Pow(radius[sections.Length - 1], 2));
            Console.WriteLine("sigma : " + sigma);

            double dy = 0.0;
            dy = (radius[sections.Length - 1] - radius[0]) / sections.Length;
            Console.WriteLine("dy : " + dy);

            double Cl_alpha = 2 * Math.PI;
            Console.WriteLine("Cl_alpha : " + Cl_alpha);

            double B = 0.95;

            double Ct = 0.0;
            Ct = Calc.Get_Ct(sigma, Cl_alpha, theta[0], B);
            Console.WriteLine("Ct : " + Ct);

            double rambda = 0.0;
            rambda = Calc.Get_rambda(Ct);
            Console.WriteLine("rambda : " + rambda);
            Console.WriteLine();

            List<double> U = new List<double>();
            List<double> phi = new List<double>();

            for (int i = 0; i < sections.Length; i++)
            {
                U.Add(Calc.Get_UT(omega, radius[i]));
                Console.WriteLine("U[" + i + "] : " + U[i]);
            }

            // "phi = rambda / r", "r = 블레이드의 상대 위치" 
            for (int i = 0; i < sections.Length; i++)
            {
                phi.Add(Calc.Get_phi(rambda, (radius[i] / radius[sections.Length - 1])));
                Console.WriteLine("phi[" + i + "] : " + phi[i]);
            }
            Console.WriteLine();


            List<double> alpha = new List<double>();
            for (int i = 0; i < sections.Length; i++)
            {
                alpha.Add(theta[i] - phi[i]);

                Console.WriteLine("alpha[" + i + "] : " + alpha[i]);
            }
            Console.WriteLine();

            ReadData Read = new ReadData();

            string file_path = @"C:\csbemt_v2\LookUpTable_NACA0012.txt";

            List<double> reynolds_LookUp = new List<double>();
            List<double> alpha_LookUp = new List<double>();
            List<List<double>> Cl_LookUp = new List<List<double>>();

            Read.Airfoil_dat(file_path, reynolds_LookUp, alpha_LookUp, Cl_LookUp);

            List<double> Cl = new List<double>();
            List<double> Cd = new List<double>();
            List<double> alpha_calc = new List<double>();
            List<double> reynolds_calc = new List<double>();

            // dat파일에서 twist에 가장 근접한 alpha 찾기 
            List<int> alpha_index = new List<int>();
            List<int> reynolds_index = new List<int>();

            int index = 0;
            for (int i = 0; i < sections.Length; i++)
            {
                double alpha_ = 0;
                for (int j = 0; j < alpha_LookUp.Count - 1; j++)
                {
                    if (Math.Abs(alpha[i] - alpha_LookUp[j]) > Math.Abs(alpha[i] - alpha_LookUp[j + 1]))
                    {
                        alpha_ = alpha_LookUp[j + 1];
                        index++;
                    }
                    else
                        break;
                }
                alpha_calc.Add(alpha_);
                alpha_index.Add(index);
            }

            index = 0;
            for (int i = 0; i < sections.Length; i++)
            {
                double reynolds_ = 0.0;
                for (int j = 0; j < reynolds_LookUp.Count; j++)
                {
                    if (Math.Abs(reynolds[i] - reynolds_LookUp[j]) > Math.Abs(reynolds[i] - reynolds_LookUp[j + 1]))
                    {
                        reynolds_ = reynolds_LookUp[j + 1];

                    }
                    else
                        break;
                }

                reynolds_calc.Add(reynolds_);
                reynolds_index.Add(index);
            }

            for (int i = 0; i < alpha_calc.Count; i++)
            {
                Console.WriteLine("alpha_calc[" + i + "] :" + alpha_calc[i]);
            }
            Console.WriteLine();

            for (int i = 0; i < reynolds_calc.Count; i++)
            {
                Console.WriteLine("reynolds_calc[" + i + "] :" + reynolds_calc[i]);
            }
            Console.WriteLine();

            Console.WriteLine(reynolds[reynolds_index[0]]);

            for (int i = 0; i < 361; i++)
            {
                List<double> test = new List<double>();
                test.Add(Calc.Interpolation(reynolds[reynolds_index[i]],
                                            reynolds[i], reynolds[i+1],
                                            Cl_LookUp[i+1][reynolds_index[i]],
                                            Cl_LookUp[i + 1][reynolds_index[i+1]]));
            }


            //double test = Calc.Interpolation(reynolds[0], reynolds_LookUp[0], reynolds_LookUp[1], );
            //Console.WriteLine(test);

            //double Lift = 0.0;
            //List<double> dL = new List<double>();
            //for (int i = 0; i < sections.Length; i++)
            //{
            //    dL.Add(Calc.Get_dL(rho, U[i], chord[i], Cl[i], dy));
            //    Console.WriteLine("dL[" + i + "] : " + dL[i]);

            //    Lift += dL[i];
            //}
            //Console.WriteLine();

            //double Drag = 0.0;
            //List<double> dD = new List<double>();
            //for (int i = 0; i < sections.Length; i++)
            //{
            //    dD.Add(Calc.Get_dD(rho, U[i], chord[i], Cd[i], dy));
            //    Console.WriteLine("dD[" + i + "] : " + dD[i]);

            //    Drag += dD[i];
            //}
            //Console.WriteLine();

            //List<double> dFx = new List<double>();
            //for (int i = 0; i < sections.Length; i++)
            //{
            //    dFx.Add(Calc.Get_dFx(dL[i], dD[i], phi[i]));
            //    Console.WriteLine("dFx[" + i + "] : " + dFx[i]);
            //}
            //Console.WriteLine();

            //List<double> dFz = new List<double>();
            //for (int i = 0; i < sections.Length; i++)
            //{
            //    dFz.Add(Calc.Get_dFz(dL[i], dD[i], phi[i]));
            //    Console.WriteLine("dFz[" + i + "] : " + dFz[i]);
            //}
            //Console.WriteLine();


            //double Thrust = 0.0;
            //List<double> dT = new List<double>();
            //for (int i = 0; i < sections.Length; i++)
            //{
            //    dT.Add(Calc.Get_dT(Nb, dFz[i]));
            //    Console.WriteLine("dT[" + i + "] : " + dT[i]);

            //    Thrust += 2 * (Math.PI) * radius[i] * dT[i];
            //}
            //Console.WriteLine();


            //double Torque = 0.0;
            //List<double> dQ = new List<double>();
            //for (int i = 0; i < sections.Length; i++)
            //{
            //    dQ.Add(Calc.Get_dQ(dT[i], radius[i]));
            //    Console.WriteLine("dQ[" + i + "] : " + dQ[i]);

            //    Torque += 2 * (Math.PI) * radius[i] * dQ[i];
            //}
            //Console.WriteLine();

            //double Power = 0.0;
            //List<double> dP = new List<double>();
            //for (int i = 0; i < sections.Length; i++)
            //{
            //    dP.Add(Calc.Get_dP(dQ[i], omega));
            //    Console.WriteLine("dP[" + i + "] : " + dP[i]);

            //    Power += 2 * (Math.PI) * radius[i] * dP[i];
            //}
            //Console.WriteLine();

            //Console.WriteLine("Lift : " + Lift);
            //Console.WriteLine("Drag : " + Drag);
            //Console.WriteLine("Thrust : " + Thrust);
            //Console.WriteLine("Torque : " + Torque);
            //Console.WriteLine("Power : " + Power);



            //double test = GaussLegendreRule.Integrate( => 9, 0.0, 2 * Math.PI, 10);
            //Console.WriteLine(test);

            //double Thrust_Total = 0;
            //List<double> Thrust = new List<double>();

            //for (int i = 0; i < swections.Length; i++)
            //{
            //    Thrust.Add(Calculation.Get_T(Ct[i], rho, radius[i], omega));

            //    Thrust_Total += Thrust[i];

            //    Console.WriteLine("T[" + i + "] : " + Math.Round(Thrust[i], 3));

            //}

            //Console.WriteLine("Thrust_Total : " + Math.Round(Thrust_Total, 3));

            //double Lift = 0.0;

            //List<double> dL = new List<double>();

            //for (int i = 0; i < sections.Length-1; i++)
            //{
            //    dL.Add(Calculation.Get_dL(rho, omega, radius[i], chord[i], Cl[i], radius[i+1]));
            //    Console.WriteLine("dL[" + i + "] : " + dL[i]);

            //    Lift += dL[i];
            //}

            //Console.WriteLine("Lift : " + Lift);
        }

    }
}
