using System;
using System.Collections.Generic;
using System.Dynamic;
using System.Linq;
using System.Text;
using System.Text.Encodings.Web;
using System.Threading.Tasks;
using static System.Net.Mime.MediaTypeNames;

namespace csbemt_v2
{
    internal class Program
    {
        const double PI = 3.14159265358979323846;
        const double deg2rad = PI / 180;
        const double rad2deg = 180 / PI;
        const double g = 9.81;
        const double rpm2omega = 2 * PI / 60;
        const double omega2rpm = 60 / (2 * PI);

        static void Main()
        {
            Calculate Calc = new Calculate();
            ReadWriteINIfile ini = new ReadWriteINIfile(@"C:\csbemt_v2\rotor.ini");

            ////////////////////
            /// ini 파일 읽기 ///
            ////////////////////
            //case
            double rpm = double.Parse(ini.ReadINI("case", "rpm"));
            double v_inf = double.Parse(ini.ReadINI("case", "v_inf"));

            //rotor
            double nblades = double.Parse(ini.ReadINI("rotor", "nblades"));
            double diameter = double.Parse(ini.ReadINI("rotor", "diameter"));
            double radius_hub = double.Parse(ini.ReadINI("rotor", "radius_hub"));
            string[] sections = ini.ReadINI("rotor", "section").Split(' ');

            ////string으로 불러오기
            string[] str_radius = ini.ReadINI("rotor", "radius").Split(' ');
            string[] str_chord = ini.ReadINI("rotor", "chord").Split(' ');
            string[] str_theta = ini.ReadINI("rotor", "pitch").Split(' ');

            ////string -> double 변환
            List<double> radius = new List<double>();
            List<double> chord = new List<double>();
            List<double> theta = new List<double>();

            for (int i = 0; i < sections.Length; i++)
            {
                radius.Add(double.Parse(str_radius[i]));
                chord.Add(double.Parse(str_chord[i]));
                theta.Add(double.Parse(str_theta[i]));
            }

            if (sections.Length != radius.Count)
                Console.WriteLine("radius.Count 오류\n");
            else if (sections.Length != chord.Count)
                Console.WriteLine("chord.Count 오류\n");
            else if (sections.Length != theta.Count)
                Console.WriteLine("theta.Count 오류\n");

            //fluid
            double rho = double.Parse(ini.ReadINI("fluid", "rho"));
            double mu = double.Parse(ini.ReadINI("fluid", "mu"));

            ////////////////////
            /// ini 정보 출력 ///
            ////////////////////
            Console.WriteLine("[case]");
            Console.WriteLine("rpm = " + rpm + "");
            Console.WriteLine("v_inf = " + v_inf + "\n");
            Console.WriteLine("[rotor]");
            Console.WriteLine("nblades = " + nblades + "");
            Console.WriteLine("diameter = " + diameter + "");
            Console.WriteLine("radius_hub = " + radius_hub + "");
            Console.WriteLine("section = " + string.Join(" ", sections) + "");
            Console.WriteLine("radius = " + string.Join(" ", radius) + "");
            Console.WriteLine("chord = " + string.Join(" ", chord) + "");
            Console.WriteLine("pitch = " + string.Join(" ", theta) + "\n");
            Console.WriteLine("[fluid]");
            Console.WriteLine("rho = " + rho +  "");
            Console.WriteLine("mu = " + mu + "\n");

            ////////////////////////////
            /// omega, reynolds 계산 ///
            ////////////////////////////
            //omega
            double omega = rpm * rpm2omega;
            Console.WriteLine("omega = " + omega.ToString("N2"));

            //reynolds
            List<double> reynolds = new List<double>();
            for (int i = 0; i < sections.Length; i++)
            {
                reynolds.Add(Calc.Get_Reynolds(rho, omega, radius[i], chord[i], mu));
                Console.WriteLine("reynolds[" + i + "] : " + reynolds[i].ToString("N2"));
            }
            //Tostring에서 숫자 쉼표 제거하는 법?



            //Console.WriteLine();
            double sigma = 0.0;
            sigma = (nblades * chord[sections.Length - 1] * (radius[sections.Length - 1] - radius[0])) / (Math.PI * Math.Pow(radius[sections.Length - 1], 2));
            //Console.WriteLine("sigma : " + sigma);

            double dy = 0.0;
            dy = (radius[sections.Length - 1] - radius[0]) / sections.Length;
            //Console.WriteLine("dy : " + dy);

            double Cl_alpha = 2 * Math.PI;
            //Console.WriteLine("Cl_alpha : " + Cl_alpha);

            double B = 0.95;

            double Ct = 0.0;
            Ct = Calc.Get_Ct(sigma, Cl_alpha, theta[0], B);
            //Console.WriteLine("Ct : " + Ct);

            double rambda = 0.0;
            rambda = Calc.Get_rambda(Ct);
            //Console.WriteLine("rambda : " + rambda);
            //Console.WriteLine();

            List<double> U = new List<double>();
            List<double> phi = new List<double>();

            for (int i = 0; i < sections.Length; i++)
            {
                U.Add(Calc.Get_UT(omega, radius[i]));
                //Console.WriteLine("U[" + i + "] : " + U[i]);
            }

            // "phi = rambda / r", "r = 블레이드의 상대 위치" 
            for (int i = 0; i < sections.Length; i++)
            {
                phi.Add(Calc.Get_phi(rambda, (radius[i] / radius[sections.Length - 1])));
                //Console.WriteLine("phi[" + i + "] : " + phi[i]);
            }
            //Console.WriteLine();


            List<double> alpha = new List<double>();
            for (int i = 0; i < sections.Length; i++)
            {
                alpha.Add(theta[i] - phi[i]);

                //Console.WriteLine("alpha[" + i + "] : " + alpha[i]);
            }
            //Console.WriteLine();

            ReadData Read = new ReadData();

            string file_path = @"C:\csbemt_v2\LookUpTable_NACA0012.txt";

            List<double> reynolds_LookUp = new List<double>();
            List<double> alpha_LookUp = new List<double>();
            List<List<double>> Cd_LookUp = new List<List<double>>();
            Read.Airfoil_dat(file_path, reynolds_LookUp, alpha_LookUp, Cd_LookUp);

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
                index = 0;
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
                        index++;

                    }
                    else
                        break;
                }

                reynolds_calc.Add(reynolds_);
                reynolds_index.Add(index);
                index = 0;
            }

            for (int i = 0; i < alpha_calc.Count; i++)
            {
                //Console.WriteLine("alpha_calc[" + i + "] :" + alpha_calc[i]);
            }
            //Console.WriteLine();

            for (int i = 0; i < reynolds_calc.Count; i++)
            {
                //Console.WriteLine("reynolds_calc[" + i + "] :" + reynolds_calc[i]);
            }
            //Console.WriteLine();

            //Console.WriteLine(reynolds[reynolds_index[0]]);

            //Console.WriteLine("alpha 배열 크기 : " + alpha.Count);                         // 10
            //Console.WriteLine("alpha_LookUp 배열 크기 : " + alpha_LookUp.Count);           // 361
            //Console.WriteLine("alpha_index 배열 크기 : " + alpha_index.Count);             // 10

            //Console.WriteLine("reynolds 배열 크기 : " + reynolds.Count);                   // 10
            //Console.WriteLine("reynolds_LookUp 배열 크기 : " + reynolds_LookUp.Count);     // 20
            //Console.WriteLine("reynolds_index 배열 크기 : " + reynolds_index.Count);       // 10

            //Console.WriteLine("Cl_LookUp 배열 세로 크기 : " + Cd_LookUp.Count);            // 362
            //Console.WriteLine("Cl_LookUp 배열 가로 크기 : " + Cd_LookUp[1].Count);         // 20

            for (int i = 0; i < alpha_index.Count; i++)
            {
                //Console.WriteLine("alpha_index[" + i + "] : " + alpha_index[i]);
            }

            for (int i = 0; i < reynolds_index.Count; i++)
            {
                //Console.WriteLine("reynolds_index[" + i + "] : " + reynolds_index[i]);
            }

            List<double> Cl = new List<double>();

            for (int i = 0; i < sections.Length; i++)
            {
                if (alpha[i] > alpha_calc[i])
                {
                    double test = Calc.Interpolation(alpha[i],
                                                     alpha_LookUp[alpha_index[i]],
                                                     alpha_LookUp[alpha_index[i] + 1],
                                                     Cd_LookUp[alpha_index[i] + 1][reynolds_index[i]],
                                                     Cd_LookUp[alpha_index[i] + 2][reynolds_index[i]]);
                    //Console.WriteLine("test " + i + ": " + test);
                }
                else
                {
                    double test1 = Calc.Interpolation(alpha[i],
                                                     alpha_LookUp[alpha_index[i]],
                                                     alpha_LookUp[alpha_index[i] - 1],
                                                     Cd_LookUp[alpha_index[i] + 1][reynolds_index[i]],
                                                     Cd_LookUp[alpha_index[i]][reynolds_index[i]]);

                    double test2 = Calc.Interpolation(alpha[i],
                                                   alpha_LookUp[alpha_index[i]],
                                                   alpha_LookUp[alpha_index[i] - 1],
                                                   Cd_LookUp[alpha_index[i] + 1][reynolds_index[i] + 1],
                                                   Cd_LookUp[alpha_index[i]][reynolds_index[i] + 1]);

                    if (reynolds[i] > reynolds_calc[i])
                    {
                        double test3 = Calc.Interpolation(reynolds[i],
                                                          reynolds_LookUp[reynolds_index[i]],
                                                          reynolds_LookUp[reynolds_index[i] + 1],
                                                          test1, test2);

                        Cl.Add(test3);
                        //Console.WriteLine("test1 " + i + ": " + test1);
                        //Console.WriteLine("test2 " + i + ": " + test2);
                        //Console.WriteLine("test3 " + i + ": " + test3);
                        //Console.WriteLine("Cl[" + i + "] :" + Cl[i]);
                        //Console.WriteLine("Cl.count : " + Cl.Count);
                        //Console.WriteLine();
                    }
                    else
                    {
                        double test3 = Calc.Interpolation(reynolds[i],
                                                          reynolds_LookUp[reynolds_index[i]],
                                                          reynolds_LookUp[reynolds_index[i] - 1],
                                                          test2, test1);

                        Cl.Add(test3);
                        //Console.WriteLine("test1 " + i + ": " + test1);
                        //Console.WriteLine("test2 " + i + ": " + test2);
                        //Console.WriteLine("test3 " + i + ": " + test3);
                        //Console.WriteLine("Cl[" + i + "] :" + Cl[i]);
                        //Console.WriteLine("Cl.count : " + Cl.Count);
                        //Console.WriteLine();
                    }
                }

                //Console.WriteLine("alpha[" + i + "] : " + alpha[i]); 
                //Console.WriteLine("alpha_LookUp[alpha_index[" + i + "]] : " + alpha_LookUp[alpha_index[i]]);
                //Console.WriteLine("alpha_LookUp[alpha_index[" + i + "]-1] : " + alpha_LookUp[alpha_index[i] - 1]);  ;
                //Console.WriteLine("Cl_LookUp[alpha_index[" + i + "] + 1][reynolds_index[" + i + "]] " + Cl_LookUp[alpha_index[i] + 1][reynolds_index[i]]);
                //Console.WriteLine("Cl_LookUp[alpha_index[" + i + "]][reynolds_index[" + i + "]] : " + Cl_LookUp[alpha_index[i]][reynolds_index[i]]);


            }

            reynolds_LookUp.Clear();
            alpha_LookUp.Clear();


            Read.Airfoil_dat(file_path, reynolds_LookUp, alpha_LookUp, Cd_LookUp);
            List<double> Cd = new List<double>();
            for (int i = 0; i < sections.Length; i++)
            {
                if (alpha[i] > alpha_calc[i])
                {
                    double test = Calc.Interpolation(alpha[i],
                                                     alpha_LookUp[alpha_index[i]],
                                                     alpha_LookUp[alpha_index[i] + 1],
                                                     Cd_LookUp[alpha_index[i] + 1][reynolds_index[i]],
                                                     Cd_LookUp[alpha_index[i] + 2][reynolds_index[i]]);
                }
                else
                {
                    double test1 = Calc.Interpolation(alpha[i],
                                                     alpha_LookUp[alpha_index[i]],
                                                     alpha_LookUp[alpha_index[i] - 1],
                                                     Cd_LookUp[alpha_index[i] + 1][reynolds_index[i]],
                                                     Cd_LookUp[alpha_index[i]][reynolds_index[i]]);

                    double test2 = Calc.Interpolation(alpha[i],
                                                   alpha_LookUp[alpha_index[i]],
                                                   alpha_LookUp[alpha_index[i] - 1],
                                                   Cd_LookUp[alpha_index[i] + 1][reynolds_index[i] + 1],
                                                   Cd_LookUp[alpha_index[i]][reynolds_index[i] + 1]);

                    if (reynolds[i] > reynolds_calc[i])
                    {
                        double test3 = Calc.Interpolation(reynolds[i],
                                                          reynolds_LookUp[reynolds_index[i]],
                                                          reynolds_LookUp[reynolds_index[i] + 1],
                                                          test1, test2);

                        Cd.Add(test3);

                        //Console.WriteLine("test1 " + i + ": " + test1);
                        //Console.WriteLine("test2 " + i + ": " + test2);
                        //Console.WriteLine("test3 " + i + ": " + test3);
                        //Console.WriteLine();
                    }
                    else
                    {
                        double test3 = Calc.Interpolation(reynolds[i],
                                                          reynolds_LookUp[reynolds_index[i]],
                                                          reynolds_LookUp[reynolds_index[i] - 1],
                                                          test2, test1);

                        Cd.Add(test3);
                        //Console.WriteLine("test1 " + i + ": " + test1);
                        //Console.WriteLine("test2 " + i + ": " + test2);
                        //Console.WriteLine("test3 " + i + ": " + test3);
                        //Console.WriteLine();
                    }
                }

                //Console.WriteLine("alpha[" + i + "] : " + alpha[i]); 
                //Console.WriteLine("alpha_LookUp[alpha_index[" + i + "]] : " + alpha_LookUp[alpha_index[i]]);
                //Console.WriteLine("alpha_LookUp[alpha_index[" + i + "]-1] : " + alpha_LookUp[alpha_index[i] - 1]);  ;
                //Console.WriteLine("Cl_LookUp[alpha_index[" + i + "] + 1][reynolds_index[" + i + "]] " + Cl_LookUp[alpha_index[i] + 1][reynolds_index[i]]);
                //Console.WriteLine("Cl_LookUp[alpha_index[" + i + "]][reynolds_index[" + i + "]] : " + Cl_LookUp[alpha_index[i]][reynolds_index[i]]);


            }

            double Lift = 0.0;
            List<double> dL = new List<double>();
            for (int i = 0; i < sections.Length; i++)
            {
                dL.Add(Calc.Get_dL(rho, U[i], chord[i], Cl[i], dy));
                //Console.WriteLine("dL[" + i + "] : " + dL[i]);

                Lift += dL[i];
            }
            //Console.WriteLine();

            double Drag = 0.0;
            List<double> dD = new List<double>();
            for (int i = 0; i < sections.Length; i++)
            {
                dD.Add(Calc.Get_dD(rho, U[i], chord[i], Cd[i], dy));
                //Console.WriteLine("dD[" + i + "] : " + dD[i]);

                Drag += dD[i];
            }
            //Console.WriteLine();

            List<double> dFx = new List<double>();
            for (int i = 0; i < sections.Length; i++)
            {
                dFx.Add(Calc.Get_dFx(dL[i], dD[i], phi[i]));
                //Console.WriteLine("dFx[" + i + "] : " + dFx[i]);
            }
            //Console.WriteLine();

            List<double> dFz = new List<double>();
            for (int i = 0; i < sections.Length; i++)
            {
                dFz.Add(Calc.Get_dFz(dL[i], dD[i], phi[i]));
                //Console.WriteLine("dFz[" + i + "] : " + dFz[i]);
            }
            //Console.WriteLine();


            double Thrust = 0.0;
            List<double> dT = new List<double>();
            for (int i = 0; i < sections.Length; i++)
            {
                dT.Add(Calc.Get_dT(nblades, dFz[i]));
                //Console.WriteLine("dT[" + i + "] : " + dT[i]);

                Thrust += 2 * (Math.PI) * radius[i] * dT[i];
            }
            //Console.WriteLine();


            double Torque = 0.0;
            List<double> dQ = new List<double>();
            for (int i = 0; i < sections.Length; i++)
            {
                dQ.Add(Calc.Get_dQ(dT[i], radius[i]));
                //Console.WriteLine("dQ[" + i + "] : " + dQ[i]);

                Torque += 2 * (Math.PI) * radius[i] * dQ[i];
            }
            //Console.WriteLine();

            double Power = 0.0;
            List<double> dP = new List<double>();
            for (int i = 0; i < sections.Length; i++)
            {
                dP.Add(Calc.Get_dP(dQ[i], omega));
                //Console.WriteLine("dP[" + i + "] : " + dP[i]);

                Power += 2 * (Math.PI) * radius[i] * dP[i];
            }
            //Console.WriteLine();

            //Console.WriteLine("Lift : " + Lift);
            //Console.WriteLine("Drag : " + Drag);
            Console.WriteLine("Thrust : " + Thrust);
            Console.WriteLine("Torque : " + Torque);
            Console.WriteLine("Power : " + Power);



        }

    }
}
