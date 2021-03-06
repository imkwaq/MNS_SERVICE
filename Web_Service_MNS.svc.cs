﻿using System;

namespace MNS_SERVICE
{
    // ПРИМЕЧАНИЕ. Команду "Переименовать" в меню "Рефакторинг" можно использовать для одновременного изменения имени класса "Web_Service_MNS" в коде, SVC-файле и файле конфигурации.
    // ПРИМЕЧАНИЕ. Чтобы запустить клиент проверки WCF для тестирования службы, выберите элементы Web_Service_MNS.svc или Web_Service_MNS.svc.cs в обозревателе решений и начните отладку.
    public class Web_Service_MNS : IWeb_Service_MNS
    {

        #region description of the data
        const int MF = 20;
        public int nv, n, nr, nc, nl, nju, ntri, nji, nei, neu, nou, ntr, noui, ntu, ntb, nf, lp, lm, kp, km, k;

        public float[] f = new float[MF + 1];
        public float[] kum = new float[MF + 1];
        public float[] kua = new float[MF + 1];
        public float[] rim = new float[MF + 1];
        public float[] ria = new float[MF + 1];
        public float[] rom = new float[MF + 1];
        public float[] roa = new float[MF + 1];

        public int[,] in_r;
        public float[] z_r;
        public int[] In_r;

        public int[,] in_c;
        public float[] z_c;
        public int[] In_c;

        public int[,] in_l;
        public float[] z_l;
        public int[] In_l;

        public int[,] in_ju;
        public int[] In_ju;
        public float[,] z_ju;

        public int[,] in_eu;
        public int[] In_eu;
        public float[,] z_eu;       

        public float[,] z_ji;
        public int[,] in_ji;
        public int[] In_ji;

        public float[,] z_ei;
        public int[,] in_ei;
        public int[] In_ei;

        public int[,] in_ou;
        public int[] In_ou;
        public float[,] z_ou;

        public int[,] in_oui;
        public int[] In_oui;
        public float[] z_oui;

        public int[,] in_tu;
        public int[] In_tu;
        public float[,] z_tu;

        public int[,] in_tb;
        public int[] In_tb;
        public float[,] z_tb;

        public int[,] in_tri;
        public int[] In_tri;
        public float[] z_tri;

        public int[,] in_tr;
        public int[] In_tr;
        public float[,] z_tr;

        public float[] Out = new float[6 * MF + 1];

        #endregion
        private int[,] get2Darray(int[] array, int nodes, int count) 
        {
            int[,] result = new int[count + 1, nodes];
            for (int i = 0; i < count; i++)
                for (int j = 0; j < nodes; j++)
                    result[i + 1, j] = array[i * nodes + j];
            return result;
        }
        private float[,] get2Darray(float[] array, int nodes, int count)
        {
            float[,] result = new float[count + 1, nodes];
            for (int i = 0; i < count; i++)
                for (int j = 0; j < nodes; j++)
                    result[i + 1, j] = array[i * nodes + j];
            return result;
        }

        public float[] OnCalc(
            int[] In_r, float[] z_r, int nr, 
            int[] In_c, float[] z_c, int nc, 
            int[] In_l, float[] z_l, int nl, 
            int nv, int lp, int lm, int kp, int km, float[] f, int nf, 
            int nju, int[] In_ju, float[] Z_ju, 
            int neu, int[] In_eu, float[] Z_eu, 
            int nji, int[] In_ji, float[] Z_ji,
            int nei, int[] In_ei, float[] Z_ei,
            int nou, int[] In_ou, float[] Z_ou, 
            int ntri, int[] In_tri, float[] Z_tri,
            int ntb, int[] In_tb, float[] Z_tb,
            int ntu, int[] In_tu, float[] Z_tu,
            int ntr, int[] In_tr, float[] Z_tr,
            int noui, int[] In_oui, float[] Z_oui)
        {
            #region Unpacking one-dimensional arrays
            in_r = get2Darray(In_r, 2, nr);
            in_c = get2Darray(In_c, 2, nc);
            in_l = get2Darray(In_l, 2, nl);

            in_ju = get2Darray(In_ju, 4, nju);
            z_ju = get2Darray(Z_ju, 3, nju);
            in_eu = get2Darray(In_eu, 4, neu);
            z_eu = get2Darray(Z_eu, 3, neu);
            in_ji = get2Darray(In_ji, 4, nji);
            z_ji = get2Darray(Z_ji, 3, nji);
            in_ei = get2Darray(In_ei, 4, nei);
            z_ei = get2Darray(Z_ei, 3, nei);
        
            in_tu = get2Darray(In_tu, 4, ntu);
            z_tu = get2Darray(Z_tu, 5, ntu);
            in_tb = get2Darray(In_tb, 6, ntb);
            z_tb = get2Darray(Z_tb, 4, ntb);
            in_tr = get2Darray(In_tr, 4, ntr);
            z_tr = get2Darray(Z_tr, 5, ntr);
            in_tri = get2Darray(In_tri, 4, ntri);
            z_tri = Z_tri;
            in_ou = get2Darray(In_ou, 5, nou);
            z_ou = get2Darray(Z_ou, 4, nou);
            in_oui = get2Darray(In_oui, 4, noui);
            z_oui = Z_oui;
            #endregion

            for (int kf = 1; kf <= nf; kf++)
            {
                n = nv + nr + nc + nl + nju + neu + nji + nou + ntri + ntr + ntu + ntb + noui; 
                Complex[,] w = new Complex[n + 1, n + 1];
                Complex cn = new Complex(0, 0);
                for (int i = 0; i <= n; i++)
                    for (int j = 0; j <= n; j++)
                    {
                        w[i, j] = cn;
                    }
                Complex s = new Complex(0.0, 2 * 3.141593 * f[kf]);
                n = nv;

                form_d(in_r, z_r, nr, 'R', w, s);
                form_d(in_c, z_c, nc, 'C', w, s);
                form_l(in_l, z_l, nl, w, s);
                
                form_ju(in_ju, z_ju, nju, w, s);
                form_ji(in_ji, z_ji, nji, w, s);
                form_eu(in_eu, z_eu, neu, w, s);
                form_ei(in_ei, z_ei, nei, w, s);

                form_tr(in_tb, z_tb, ntb, w, s);
                form_tu(in_tu, z_tu, ntu, w, s);
                form_tr(in_tr, z_tr, ntr, w, s);
                form_tri(in_tri, z_tri, ntri, w);                
                form_ou(in_ou, z_ou, nou, w, s);
                form_oui(in_oui, z_oui, noui, w);

                form_s(lp, lm, w);

                if ((lp == 1) && (lm == 0) && (kp == 2) && (km == 0))
                {
                    st(w, n);
                    sf1(kf, w, ref kum, ref kua, ref rim, ref ria, ref rom, ref roa);
                }
                else
                {
                    gauss_c(w, n);
                    sf2(kf, w, lp, lm, kp, km, ref kum, ref kua, ref rim, ref ria);
                }
            }

            #region Packing Output Array - Out
            for (int i = 1; i <= MF; i++)
            {
                Out[i] = kum[i];
                Out[i + MF] = kua[i];
                Out[i + 2 * MF] = rim[i];
                Out[i + 3 * MF] = ria[i];
                Out[i + 4 * MF] = rom[i];
                Out[i + 5 * MF] = roa[i];
            }
#endregion
            return Out;
        }

        //Формирование комплексных частных матриц двхполюсников R, L, C
        public void form_d(int[,] in_d, float[] z_d, int nd, char td, Complex[,] w, Complex s)
        {
            int i, j, g;
            for (int kd = 1; kd <= nd; kd++)
                for (int l = 0; l <= 1; l++)
                {
                    i = in_d[kd, l];
                    if (i == 0) continue;
                    for (int m = 0; m <= 1; m++)
                    {
                        j = in_d[kd, m];
                        if (j == 0) continue;
                        g = (1 - 2 * l) * (1 - 2 * m);
                        switch (td)
                        {
                            case 'R': w[i, j] += g / z_d[kd]; break;
                            case 'C': w[i, j] += g * s * z_d[kd]; break;
                                //case 'L': w[i, j] += g / (s * z_d[kd]); break;
                        }
                    }
                }
        }

        public void form_l(int[,] in_d, float[] z_d, int nd, Complex[,] w, Complex s)
        {
            int i, j, g;
            for (int kd = 1; kd <= nd; kd++)
            {
                i = n + kd;
                w[i, i] = s * z_d[kd];
                for (int m = 0; m <= 1; m++)
                {
                    j = in_d[kd, m];
                    if (j == 0) continue;
                    g = (1 - 2 * m);
                    w[i, j] -= g;
                    w[j, i] += g;
                }
            }
            n += nd;
        }

        //Формирование комплексных частных матриц частотно-зависимых ИТУН
        public void form_ju(int[,] in_ju, float[,] z_ju, int nju, Complex[,] w, Complex s)
        {
            Complex ys = new Complex(0, 0);
            int i, j, g;
            for (int kju = 1; kju <= nju; kju++)
            {
                ys = z_ju[kju, 0] * (1 + s * z_ju[kju, 1]) / (1 + s * z_ju[kju, 2]);
                for (int l = 2; l <= 3; l++)
                {
                    i = in_ju[kju, l];
                    if (i == 0) continue;
                    for (int m = 0; m <= 1; m++)
                    {
                        j = in_ju[kju, m];
                        if (j == 0) continue;
                        g = (5 - 2 * l) * (1 - 2 * m);
                        w[i, j] += g * ys;
                    }
                }
            }
        }        

        //Формирование конмлексных частных матриц частотно-зависимого ИНУН
        public void form_eu(int[,] in_eu, float[,] z_eu, int neu, Complex[,] w, Complex s)
        {
            Complex nu = new Complex(0, 0);
            int i, j, g;
            for (int keu = 1; keu <= neu; keu++)
            {
                nu = z_eu[keu, 0] * (1+ s * z_eu[keu, 1]) / (1+ s * z_eu[keu, 2]);
                i = n + keu;
                for (int m = 0; m <= 3; m++)
                {
                    j = in_eu[keu, m];
                    if (j == 0) continue;
                    if (m < 2)
                    {
                        g = 1 - 2 * m;
                        w[i, j] += g * nu;
                    }
                    else
                    {
                        g = 5 - 2 * m;
                        w[i, j] -= g;
                        w[j, i] += g;
                    }
                }
            }
            n += neu;
        }

        //Формирование конмлексных частных матриц частотно-зависимого ИТУТ
        public void form_ji(int[,] in_ji, float[,] z_ji, int nji, Complex[,] w, Complex s)
        {
            Complex b = new Complex(0, 0);
            int i, j, g;
            for (int kji = 1; kji <= nji; kji++)    //цикл по всем строкам
            {
                b = z_ji[kji, 0] * (1 + s * z_ji[kji, 1]) / (1 + s * z_ji[kji, 2]);
                //for (int l = 2; l <= 3; l++)
                //{
                j = n + kji;
                for (int l = 0; l <= 3; l++)
                {
                    i = in_ji[kji, l];
                    if (i == 0) continue;
                    if (l < 2)
                    {
                        g = 1 - 2 * l;
                        w[i,j] += g;
                        w[j,i] -= g;
                    }
                    else
                    {
                        g = 5 - 2 * l;
                        w[i,j] += g * b;
                    }                    
                }
            }
        }
        //Формирование конмлексных частных матриц частотно-зависимого ИНУТ
        public void form_ei(int[,] in_ei, float[,] z_ei, int nei, Complex[,] w, Complex s)
        {
            Complex yy = new Complex(0, 0);
            int i1, i2, j, g;
            for (int kei = 1; kei <= nei; kei++)    //цикл по всем строкам
            {
                yy = z_ei[kei, 0] * (1 + s * z_ei[kei, 1]) / (1 + s * z_ei[kei, 2]);
                //for (int l = 2; l <= 3; l++)
                //{
                i1 = n + kei;
                i2 = i1 + nei;
                w[i2, i1] = yy;
                for (int m = 0; m <= 3; m++)
                {
                    j = in_ei[kei, m];
                    if (j == 0) continue;
                    if (m < 2)
                    {
                        g = 1 - 2 * m;
                        w[i1, j] -= g;
                        w[j, i1] += g;
                    }
                    else
                    {
                        g = 5 - 2 * m;
                        w[i2, j] -= g;
                        w[j, i2] += g;
                    }
                }
            }
            n += 2 * nei;
        }
        //Формирование конмлексных частных матриц идеального операционного усилителя
        public void form_oui(int[,] in_oui, float[] z_oui, int noui, Complex[,] w)
        {
            Complex yy = new Complex(0, 0);
            int i, j1, j2, g;
            for (int koui = 1; koui <= noui; koui++)    //цикл по всем строкам
            {
                j1 = n + koui;
                j2 = j1 + noui;
                w[j2, j1] = new Complex(-1, 0);
                for (int l = 0; l <= 3; l++)
                {
                    i = in_oui[koui, l];
                    if (i == 0) continue;
                    if (l < 2)
                    {
                        g = 1 - 2 * l;
                        w[i, j1] += g;
                        w[j1, i] -= g;                       
                    }
                    else
                    {
                        g = 5 - 2 * l;
                        w[i, j2] += g;
                    }
                }
            }
            n += 2 * noui;
        }
        //Формирование комплексных частных матриц трансформатора
        public void form_tr(int[,] in_tr, float[,] z_tr, int ntr, Complex[,] w, Complex s)
        {
            int i1, i2, j, g;
            for (int ktr = 1; ktr <= ntr; ktr++)
            {
                i1 = n + ktr;
                i2 = i1 + ntr;
                w[i1, i1] = z_tr[ktr, 0] + s * z_tr[ktr, 2];
                w[i2, i2] = z_tr[ktr, 1] + s * z_tr[ktr, 3];
                w[i1, i2] = s * z_tr[ktr, 4];
                w[i2, i1] = s * z_tr[ktr, 4];
                for (int m = 0; m <= 3; m++)
                {
                    j = in_tr[ktr, m];
                    if (j == 0) continue;
                    if (m < 2)
                    {
                        g = 1 - 2 * m;
                        w[i1, j] -= g;
                        w[j, i1] += g;
                    }
                    else
                    {
                        g = 5 - 2 * m;
                        w[i2, j] -= g;
                        w[j, i2] += g;
                    }
                }
            }
            n += 2 * ntr;
        }

        //Формирование комплексных частных матриц биполярного транзистора
        public void form_tb(int[,] in_tb, float[,] z_tb, int ntb, Complex[,] w, Complex s)
        {
            Complex[,] y = new Complex[5, 5];
            int i, j, ii, jj, g, l, m;
            int[,] in_d = { { 4, 3 }, { 1, 4 }, { 4, 2 }, { 1, 4 }, { 4, 2 } };
            int[] in_ju = { 4, 1, 2, 4 };
            for (int ktb = 1; ktb <= ntb; ktb++)
            {
                for (i = 1; i <= 4; i++)
                    for (j = 1; j <= 4; j++)
                        y[i, j] = new Complex();
                for (int k = 0; k <= 4; k++)
                    for (l = 0; l <= 1; l++)
                    {
                        i = in_d[k, l];
                        for (m = 0; m <= 1; m++)
                        {
                            j = in_d[k, m];
                            g = (1 - 2 * l) * (1 - 2 * m);
                            if (k < 3)
                                y[i, j] += g / z_tb[ktb, k];
                            else
                                y[i, j] += s * z_tb[ktb, k] * g;
                        }
                    }
                for (l = 2; l <= 3; l++)
                {
                    i = in_ju[l];
                    for (m = 0; m <= 1; m++)
                    {
                        j = in_ju[m];
                        g = (5 - 2 * l) * (1 - 2 * m);
                        y[i, j] += g * z_tb[ktb, 5] / (1 + z_tb[ktb, 5]) / z_tb[ktb, 1];
                    }
                }
                for (i = 3; i >= 1; i--)
                    for (j = 3; j >= 1; j--)
                        y[i, j] -= y[i, 4] * y[4, j] / y[4, 4];
                for (i = 1; i <= 3; i++)
                {
                    ii = in_tb[ktb, i];
                    if (ii == 0)
                        continue;
                    for (j = 1; j <= 3; j++)
                    {
                        jj = in_tb[ktb, j];
                        if (jj == 0)
                            continue;
                        w[ii, jj] += y[i, j];
                    }
                }                
            }
        }

        //Формирование комплексных частных матриц униполярного транзистора
        public void form_tu(int[,] in_tu, float[,] z_tu, int ntu, Complex[,] w, Complex s)
        {
            Complex[,] y = new Complex[4, 4];
            int i, j, ii, jj, g, l, m;
            int[,] in_d = { { 2, 3 }, { 1, 3 }, { 1, 2 }, { 2, 3 } };
            int[] in_ju = { 1, 3, 2, 3 };
            for (int ktu = 1; ktu <= ntu; ktu++)
            {
                for (i = 1; i <= 3; i++)
                    for (j = 1; j <= 3; j++)
                        y[i, j] = new Complex();
                for (int k = 0; k <= 3; k++)
                    for (l = 0; l <= 1; l++)
                    {
                        i = in_d[k, l];
                        for (m = 0; m <= 1; m++)
                        {
                            j = in_d[k, m];
                            g = (1 - 2 * l) * (1 - 2 * m);
                            if (k == 0)
                                y[i, j] += g / z_tu[ktu, k];
                            else
                                y[i, j] += g * s * z_tu[ktu, k];
                        }
                    }
                for (l = 2; l <= 3; l++)
                {
                    i = in_ju[l];
                    for (m = 0; m <= 1; m++)
                    {
                        j = in_ju[m];
                        g = (5 - 2 * l) * (1 - 2 * m);
                        y[i, j] += g * z_tu[ktu, 4];
                    }
                }
                
                for (i = 1; i <= 3; i++)
                {
                    ii = in_tu[ktu, i];
                    if (ii == 0)
                        continue;
                    for (j = 1; j <= 3; j++)
                    {
                        jj = in_tu[ktu, j];
                        if (jj == 0)
                            continue;
                        w[ii, jj] += y[i, j];
                    }
                }
            }
        }

        //Формирование комплексных частных матриц операционного усилителя
        public void form_ou(int[,] in_ou, float[,] z_ou, int nou, Complex[,] w, Complex s)
        {
            Complex cn = new Complex(0, 0);
            Complex[,] y = new Complex[5, 5];
            Complex ys = new Complex(0, 0);
            int[,] in_d = new int[2, 2] { { 1, 2 }, { 3, 4 } };
            int[] in_ju = new int[4] { 1, 2, 4, 3 };
            int i, j, g, ii, jj, l, m;
            for (int kou = 1; kou <= nou; kou++)
            {
                for (i = 1; i <= 4; i++)
                    for (j = 1; j <= 4; j++)
                        y[i, j] = cn;
                for (int k = 0; k <= 1; k++)
                    for (l = 0; l <= 1; l++)
                    {
                        i = in_d[k, l];
                        for (m = 0; m <= 1; m++)
                        {
                            j = in_d[k, m];
                            g = (1 - 2 * l) * (1 - 2 * m);
                            y[i, j] += g / z_ou[kou, k];
                        }
                    }
                ys = z_ou[kou, 2] / (1 + s * 0.16 * z_ou[kou, 2] /
                               z_ou[kou, 3]) / z_ou[kou, 1];
                for (l = 2; l <= 3; l++)
                {
                    i = in_ju[l];
                    for (m = 0; m <= 1; m++)
                    {
                        j = in_ju[m];
                        g = (1 - 2 * m) * (5 - 2 * l);
                        y[i, j] += g * ys;
                    }
                }
                for (i = 1; i <= 4; i++)
                {
                    ii = in_ou[kou, i];
                    if (ii == 0) continue;
                    for (j = 1; j <= 4; j++)
                    {
                        jj = in_ou[kou, j];
                        if (jj == 0) continue;
                        w[ii, jj] += y[i, j];
                    }
                }
            }
        }

        //Формирование комплексных частных матриц идеального трансформатора
        public void form_tri(int[,] in_tri, float[] z_tri, int ntri, Complex[,] w)
        {
            int i, j, g;
            for (int ktri = 1; ktri <= ntri; ktri++)
            {
                i = n + ktri;
                for (int m = 0; m <= 3; m++)
                {
                    j = in_tri[ktri, m];
                    if (j == 0) continue;
                    if (m < 2)
                    {
                        g = 1 - 2 * m;
                        w[i, j] += g * z_tri[ktri];
                        w[j, i] -= g * z_tri[ktri];
                    }
                    else
                    {
                        g = 5 - 2 * m;
                        w[i, j] -= g;
                        w[j, i] += g;
                    }
                }
            }
            n += ntri;
        }

        public void form_s(int lp, int lm, Complex[,] w)
        {
            for (int i = 0; i <= n; i++)
                w[i, 0] = new Complex(0, 0);
            if (lp != 0)
                w[lp, 0] = new Complex(-1, 0);
            if (lm != 0)
                w[lm, 0] = new Complex(1, 0);
        }

        //Приведение схемы к трехполюснику
        public void st(Complex[,] w, int n)
        {
            Complex c = new Complex(0, 0);
            Complex t = new Complex(0, 0);
            Complex cn = new Complex(0, 0);
            double g;
            int l;
            for (int k = n; k >= 3; k--)
            {
                l = k;
                g = 0.001;
                while (w[l, k].abs <= g)
                {
                    l = l - 1;
                    if (l == 2)
                    {
                        l = k;
                        g = 0.1 * g;
                    }
                }
                if (l != k)
                    for (int j = k; j <= n; j++)
                    {
                        t = w[k, j];
                        w[k, j] = w[l, j];
                        w[l, j] = t;
                    }
                for (int i = k - 1; i >= 1; i--)
                {
                    if (w[i, k] == cn)
                        continue;
                    c = w[i, k] / w[k, k];
                    for (int j = 1; j <= k - 1; j++)
                        if (w[k, j] != cn)
                            w[i, j] -= c * w[k, j];
                }
            }
        }

        //Расчет передаточных функций
        public void sf1(int kf, Complex[,] w, ref float[] kum, ref float[] kua, ref float[] rim, ref float[] ria, ref float[] rom, ref float[] roa)
        {
            Complex ku = new Complex(0, 0);
            Complex ri = new Complex(0, 0);
            Complex ro = new Complex(0, 0);
            Complex d = new Complex(0, 0);
            ku = -w[2, 1] / w[2, 2];
            d = w[1, 1] * w[2, 2] - w[2, 1] * w[1, 2];
            ri = w[2, 2] / d;
            ro = w[1, 1] / d;
            kum[kf] = (float)ku.abs;
            kua[kf] = (float)ku.arg * 180.0f / (float)Math.PI;
            rim[kf] = (float)ri.abs;
            ria[kf] = (float)ri.arg * 180.0f / (float)Math.PI;
            rom[kf] = (float)ro.abs;
            roa[kf] = (float)ro.arg * 180.0f / (float)Math.PI;
        }

        // Обобщенный метод расчета переменных
        public void gauss_c(Complex[,] w, int n)
        {
            int i, j, k, l;
            //Прямой ход
            Complex c = new Complex(0, 0);
            Complex d = new Complex(0, 0);
            Complex t = new Complex(0, 0);
            Complex cn = new Complex(0, 0);
            for (k = 1; k < n; k++)
            {
                //Перестановка строк по максимальному значению
                l = k;
                for (i = k + 1; i <= n; i++)
                    if (w[i, k].abs > w[l, k].abs)
                        l = i;
                if (l != k)
                    for (j = 0; j <= n; j++)
                        if (j == 0 || j >= k)
                        {
                            t = w[k, j];
                            w[k, j] = w[l, j];
                            w[l, j] = t;
                        }
                //Конец перестановки строк
                d = 1.0 / w[k, k];
                for (i = k + 1; i <= n; i++)
                {
                    if (w[i, k] == cn)
                        continue;
                    c = w[i, k] * d;
                    for (j = k + 1; j <= n; j++)
                        if (w[k, j] != cn)
                            w[i, j] = w[i, j] - c * w[k, j];
                    if (w[k, 0] != cn)
                        w[i, 0] = w[i, 0] - c * w[k, 0];
                }
            }
            //Обратный ход
            w[0, n] = -w[n, 0] / w[n, n];
            for (i = n - 1; i >= 1; i--)
            {
                t = w[i, 0];
                for (j = i + 1; j <= n; j++)
                    t = t + w[i, j] * w[0, j];
                w[0, i] = -t / w[i, i];
            }
        }

        //Расчет передаточных функций
        public void sf2(int kf, Complex[,] w, int lp, int lm, int kp, int km, ref float[] kum, ref float[] kua, ref float[] rim, ref float[] ria)
        {
            Complex ku = new Complex(0, 0);
            Complex ri = new Complex(0, 0);
            ri = w[0, lp] - w[0, lm];
            ku = (w[0, kp] - w[0, km]) / ri;
            kum[kf] = (float)ku.abs;
            kua[kf] = (float)ku.arg * 180.0f / (float)Math.PI;
            rim[kf] = (float)ri.abs;
            ria[kf] = (float)ri.arg * 180.0f / (float)Math.PI;
        }            
        
    }
}
