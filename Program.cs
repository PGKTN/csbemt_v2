using System;
using System.Collections.Generic;
using System.Dynamic;
using System.Globalization;
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

        internal static double Pow(double value, double exponent)
        {
            return Math.Pow(value, exponent);
        }

        internal static double Abs(double value)
        {
            return Math.Abs(value);
        }

        internal static double Sqrt(double value)
        {
            return Math.Sqrt(value);
        }

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
            Console.WriteLine("omega = " + omega.ToString("N4") + " rad/s" + "\n");

            //reynolds
            List<double> reynolds = new List<double>();
            for (int i = 0; i < sections.Length; i++)
            {
                reynolds.Add(Calc.Get_Reynolds(rho, omega, radius[i], chord[i], mu));
                Console.WriteLine("reynolds[" + i + "] : " + reynolds[i].ToString("N0").Replace(",", ""));
            }
            Console.WriteLine();

            //////////////////////
            /// dy, sigma 계산 ///
            //////////////////////
            //dy
            List<double> dy = new List<double>();
            for (int i = 0; i < sections.Length-1; i++)
            {
                dy.Add(Abs(radius[i] - radius[i+1]));
                Console.WriteLine("dy[" + i + "] = " + dy[i].ToString("N4").Replace(",", "") + " m");
            }
            Console.WriteLine();

            //A_blade, A_disk, sigma
            double A_blade = 0.0;
            for (int i = 0; i < sections.Length - 1; i++)
            {
                A_blade += chord[i] * dy[i];
            }

            double A_disk = PI * Pow(diameter / 2.0, 2);

            double sigma = A_blade / A_disk;
            Console.WriteLine("sigma = " + sigma.ToString("N4").Replace(",", "") + "\n");

            /////////////////////
            /// Cl_alpha 계산 ///
            /////////////////////
            double Cl_alpha = 2 * PI;
            Console.WriteLine("Cl_alpha = " + Cl_alpha.ToString("N4").Replace(",", "") + "\n");

            /////////////////////////
            /// Ct, U, phi, alpha ///
            /// /////////////////////
            List<double> UT = new List<double>();
            List<double> UP = new List<double>();
            List<double> U = new List<double>();
            List<double> phi = new List<double>(); 
            List<double> alpha = new List<double>();

            double B = 0.95;
            Console.WriteLine("B = " + B.ToString("N4").Replace(",", "") + "\n");

            double Ct = Calc.Get_Ct(sigma, Cl_alpha, theta[0], B);
            Console.WriteLine("Ct = " + Ct.ToString("N4").Replace(",", "") + "\n");

            double rambda = Sqrt(Ct/2.0);
            Console.WriteLine("rambda = " + rambda.ToString("N4").Replace(",", "") + "\n");

            for (int i = 0; i < sections.Length; i++)
            {
                phi.Add(rambda / radius[i]);
                UT.Add(Calc.Get_UT(omega, radius[i]));
                UP.Add(Calc.Get_UP(phi[i], UT[i]));
                U.Add(Calc.Get_U(UT[i], UP[i]));
                alpha.Add(theta[i] - phi[i]);
                Console.WriteLine("alpha[" + i + "] = " + alpha[i].ToString("N4").Replace(",", "") + " ˚");
            }
            Console.WriteLine();
        }

    }
}
