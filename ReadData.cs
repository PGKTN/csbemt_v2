using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace csbemt_v2
{
    public class ReadData
    {
        public void Airfoil_dat(string file_path, List<double> Reynolds_LookUp)
        {

            List<string> airfoil_data_line = new List<string>();

            try
            {
                // FileStream을 사용하여 파일 열기
                using (FileStream fs = new FileStream(file_path, FileMode.Open))
                {
                    // StreamReader를 사용하여 파일 읽기
                    using (StreamReader sr = new StreamReader(fs))
                    {
                        string line;

                        while ((line = sr.ReadLine()) != null)
                            airfoil_data_line.Add(line);
                    }
                }
            }
            catch (Exception ex)
            {
                // 오류 처리
                Console.WriteLine("에어포일 파일을 읽는 도중 오류가 발생했습니다: " + ex.Message);
            }

            

            for (int i = 0; i < airfoil_data_line.Count; i++)
            {
                airfoil_data_line[i] = ConsolidateTabs(airfoil_data_line[i]);

                string[] data = airfoil_data_line[i].Split('\t');

                for (int j = 0; j < data.Count; j++)
                {
                    Reynolds_LookUp.Add(data[i]);

                    Console.WriteLine(Reynolds_LookUp[j]);
                }
            }
        }

        static string ConsolidateSpaces(string input)
        {
            // 공백을 기준으로 문자열 분할
            string[] parts = input.Split(' ', StringSplitOptions.RemoveEmptyEntries);
            // 분할된 부분을 다시 합침
            return string.Join(" ", parts);
        }

        static string ConsolidateTabs(string input)
        {
            // 탭을 기준으로 문자열 분할
            string[] parts = input.Split('\t', StringSplitOptions.RemoveEmptyEntries);
            // 분할된 부분을 다시 합침
            return string.Join(" ", parts);
        }
    }
}
