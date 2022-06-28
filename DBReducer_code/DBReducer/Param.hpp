#ifndef Param_NAME
#define Param_NAME

#include <string>
#include <vector>

class Params {
public:
    virtual ~Params(){};
    Params(){};
    int tTagLen; // tag长度
    std::string selectmod; 
    std::string fixmod;
    bool isMaxQuant; // 衔接的搜索引擎是否为Maxquant
    std::string ProteinDatabase; // 原始数据库路径
    std::string SearchRes; // pFind输出结果路径
    std::string RefinedPath; // pFind输出结果路径
    double threshold = 0.01; // 谱图/肽段层次的FDR限制
    int multip; // 进程数
    std::string activation_type; // 高低精度数据区分
    int msmsnum; // 数据集数量
    int msmsfolder; // 文件夹数量
    std::vector<std::string> msmspaths; // 存储数据集的路径
    std::string msmstype; // 数据集文件类型
    void LoadParam(const std::string cfgpath);
    void PrintParam();
};


class PFindParams {
public:
    virtual ~PFindParams(){};
    PFindParams(const Params &pAnno_param) {
        multip = pAnno_param.multip;
        activation_type = pAnno_param.activation_type;
        selectmod = pAnno_param.selectmod;
        fixmod = pAnno_param.fixmod;
        psm_fdr = 0.01;
        pro_fdr = 0.01;
        open = 5; // 5是DBReducer对应的搜索模式
        modpath = pAnno_param.SearchRes + "\\modification.ini";
        fastapath = pAnno_param.ProteinDatabase;
        outputpath = pAnno_param.SearchRes;
        msmsnum = pAnno_param.msmsnum;
        msmsfolder = pAnno_param.msmsfolder;
        msmspaths.assign(pAnno_param.msmspaths.begin(), pAnno_param.msmspaths.end());
        msmstype = pAnno_param.msmstype;
    };
    PFindParams(){};
    // [param]
    int multip = 1; // 需指定
    int thread = 2;
    std::string activation_type = "HCD-FTMS"; // 需指定
    int mstol = 20;
    int mstolppm = 1;
    int msmstol = 20;
    int msmstolppm = 1;
    int temppepnum = 100;
    int pepnum = 10;
    int selectpeak = 200;
    int64_t maxprolen = 10000000000;
    int maxspec = 10000;
    int IeqL = 1;
    int npep = 2;
    int maxdelta = 500;
    std::string selectmod = ""; // 需指定
    std::string fixmod = ""; // 需指定
    int maxmod = 3;
    std::string enzyme = "Trypsin KR _ C";
    int digest = 3;
    int max_clv_sites = 3;
    // [filter]
    double psm_fdr = 0.01; // 需指定
    int psm_fdr_type = 1;
    int mass_lower = 600;
    int mass_upper = 10000;
    int len_lower = 6;
    int len_upper = 100;
    int pep_per_pro = 1;
    double pro_fdr = 0.01; // 需指定
    // [engine]
    int open = 5; // 需指定
    int open_tag_len = 5;
    int rest_tag_iteration = 1;
    int rest_tag_len = 4;
    int rest_mod_num = 10;
    int salvo_iteration = 1;
    int salvo_mod_num = 5;
    // [file]
    std::string modpath = ""; // 需指定
    std::string fastapath = ""; // 需指定
    int auto_rev = 1;
    std::string outputpath = ""; // 需指定
    // [datalist]
    int msmsnum = 0; // 需指定
    int msmsfolder = 0; // 需指定
    std::vector<std::string> msmspaths; // 需指定
    std::string msmstype = ""; // 需指定
    // [quant]
    std::string quant = "1|None";
    // [system]
    std::string log = "LOG_INFO";
    void Write(const std::string &out_path);
};
#endif