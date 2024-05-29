using System;
using System.Collections.Generic;
using System.Dynamic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace csbemt_v2
{
    public class Calculate
    {
        public double Interpolation(double x, double x1, double x2, double y1, double y2)
        {
            double y = y1 + ((x - x1) * (y2 - y1) / (x2 - x1));

            return y;
        }

        public double Bilinear_Interpolation()
        {
            return 0;
        }

        public double Get_sigma(int Nb, double chord, double radius, double root)
        {
            return (Nb * chord) * (radius) / (Math.PI * Math.Pow(radius, 2));
        }

        public double Get_Cl_alpha(double Cl_1, double Cl_2, double alpha_1, double alpha_2)
        {
            //return (Cl_1 - Cl_2) / To_Radian(alpha_1 - alpha_2);
            return 2 * Math.PI;
        }

        internal double To_Radian(double num)
        {
            return num * Math.PI / 180;
        }

        public double Get_Reynolds(double rho, double omega, double radius, double chord, double mu)
        {
            return rho * omega * radius * chord / mu;
        }


        public double Get_Ct(double sigma, double Cl_alpha, double theta, double B)
        {
            double Ct = 0;

            // 반복 식을 진행하기 위한 초기값.
            double Ct_ = 0.005;

            while (1 > 0)
            {
                Ct = 0.5 * sigma * Cl_alpha * Math.Pow(B, 2) * ((To_Radian(theta) * B / 3) - (0.5 * Math.Sqrt(Ct_ / 2)));

                if (Math.Abs(Ct - Ct_) <= Math.Pow(10, -12))
                    break;
                else
                    Ct_ = Ct;
            }

            return Ct;
        }

        public double Get_UT(double omega, double radius)
        {
            return omega * radius;
        }
        public double Get_UP(double phi, double UT)
        {
            return phi * UT;
        }

        public double Get_U(double UT, double UP)
        {
            return Math.Sqrt(Math.Pow(UT, 2) + Math.Pow(UP, 2));
        }


        //public double Get_phi(double UT, double UP)
        //{
        //    return Math.Atan(UP / UT);
        //}

        public double Get_phi(double rambda, double r)
        {
            return rambda / r;
        }

        public double Get_T(double Ct, double rho, double radius, double omega)
        {
            return Ct * rho * Math.PI * Math.Pow(radius, 2) * Math.Pow(omega * radius, 2);
        }

        public double Get_dL(double rho, double U, double chord, double Cl, double dy)
        {
            return 0.5 * rho * Math.Pow(U, 2) * chord * Cl * dy;
        }
        public double Get_dD(double rho, double U, double chord, double Cd, double dy)
        {
            return 0.5 * rho * Math.Pow(U, 2) * chord * Cd * dy;
        }

        public double Get_dFx(double dL, double dD, double phi)
        {
            return dL * Math.Sin(phi * Math.PI/180) + dD * Math.Cos(phi * Math.PI / 180);
        }

        public double Get_dFz(double dL, double dD, double phi)
        {
            return dL * Math.Cos(phi * Math.PI / 180) - dD * Math.Sin(phi * Math.PI / 180);
        }

        public double Get_rambda(double Ct)
        {
            return Math.Sqrt(Ct / 2);
        }

        public double Get_dT(double Nb, double dFz)
        {
            return Nb * dFz;
        }

        public double Get_dQ(double dT, double radius)
        {
            return dT * radius;
        }

        public double Get_dP(double DQ, double omega)
        {
            return DQ * omega;
        }
    }
}