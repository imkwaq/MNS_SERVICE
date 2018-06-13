using System.ServiceModel;

namespace MNS_SERVICE
{
    // ПРИМЕЧАНИЕ. Команду "Переименовать" в меню "Рефакторинг" можно использовать для одновременного изменения имени интерфейса "IWeb_Service_MNS" в коде и файле конфигурации.
    [ServiceContract]
    public interface IWeb_Service_MNS
    {
        [OperationContract]
        float[] OnCalc(
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
            int noui, int[] In_oui, float[] Z_oui);
    }
}
