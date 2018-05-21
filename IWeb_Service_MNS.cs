using System.ServiceModel;

namespace MNS_SERVICE
{
    // ПРИМЕЧАНИЕ. Команду "Переименовать" в меню "Рефакторинг" можно использовать для одновременного изменения имени интерфейса "IWeb_Service_MNS" в коде и файле конфигурации.
    [ServiceContract]
    public interface IWeb_Service_MNS
    {
        [OperationContract]
        float[] OnCalc(int[] In_r, float[] z_r, int nr, int[] In_c, float[] z_c, int nc, int[] In_l, float[] z_l, int nl, int nv, int lp, int lm, int kp, int km, float[] f, int nf, int nju, int[] In_ju, float[] Z_ju, int neu, int[] In_eu, float[] Z_eu, int nji, int[] In_ji, float[] Z_ji, int nou, int[] In_ou, float[] Z_ou, int ntri, int[] In_tri, float[] z_tri);
                 
            //, int njivi, int[] In_jivi, float[] z_jivi, int nju, int[] In_ju, float[] Z_ju, int ntri, int[] In_tri, float[] z_tri, int nou, int[] In_ou, float[] Z_ou
        //, GV.neu, GV.In_eu, GV.z_eu
        //int neu, int[] In_eu, float[] z_eu
    }
}
