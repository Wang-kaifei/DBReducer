#include "Param.hpp"
#include <typeinfo>
#include <iostream>
#include "../utils/IOUtil.h"
#include "../utils/StringUtil.h"
#define VNAME(value) (#value)

/**
 * @brief 从pAnno.cfg中读取参数
 * 
 * @param cfgpath pAnno.cfg文件路径
 */
void Params::LoadParam(const std::string cfgpath) {
    std::ifstream buffer = ReadFilebyStream(cfgpath);
    std::string linestr;
    while (!buffer.eof()) {
        getline(buffer, linestr);
        Strip(&linestr);
        if (linestr[0] == '#' || linestr == "" || linestr[0] == '[')
			continue;
        std::vector<std::string> segs = SplitByChars(linestr, "=#");
		std::string key = segs[0];
		std::string value = "";
		if (segs.size() > 1) {
            value = Strip(&segs[1]);
        }
		if (key == "tTagLen") 
			this->tTagLen = std::stoi(value);
        else if (key == "selectmod")
            this->selectmod = value;
        else if (key == "fixmod")
			this->fixmod = value;
        else if (key == "isMaxQuant")
			this->isMaxQuant = std::stoi(value) == 1 ? true : false;
		else if (key == "ProteinDatabase")
			this->ProteinDatabase = value;
        else if (key == "SearchRes") {
            this->SearchRes = Strip(&value, "\\/");
        }
        else if (key == "threshold")
			this->threshold = std::stod(value);
        else if (key == "multip")
			this->multip = std::stoi(value);
        else if (key == "thread")
			this->thread = std::stoi(value);
        else if (key == "maxprolen")
			this->maxprolen = std::stoll(value);
        else if (key == "maxspec")
			this->maxspec = std::stoll(value);
        else if (key == "activation_type")
			this->activation_type = value;
        else if (key == "msmsnum")
			this->msmsnum = std::stoi(value);
        else if (key == "msmsfolder")
			this->msmsfolder = std::stoi(value);
        else if (key.find("msmspath") != std::string::npos)
            this->msmspaths.emplace_back(value);
        else if (key == "msmstype")
			this->msmstype = value;
		else
			printf("Invalid param!  : %s\n", linestr.c_str());
    }
}

/**
 * @brief 临时测试用，输出读取的参数
 * 
 */
void Params::PrintParam() {
    std::cout << "tTagLen" << " = " << tTagLen << std::endl;
    std::cout << "selectmod" << " = " << selectmod << std::endl;
    std::cout << "fixmod" << " = " << fixmod << std::endl;
    std::cout << "ProteinDatabase" << " = " << ProteinDatabase << std::endl;
    std::cout << "SearchRes" << " = " << SearchRes << std::endl;
    std::cout << "RefinedPath" << " = " << RefinedPath << std::endl;
    std::cout << "threshold" << " = " << threshold << std::endl;
    std::cout << "multip" << " = " << multip << std::endl;
    std::cout << "activation_type" << " = " << activation_type << std::endl;
    std::cout << "msmsnum" << " = " << msmsnum << std::endl;
    std::cout << "msmsfolder" << " = " << msmsfolder << std::endl;
    std::cout << "msmstype" << " = " << msmstype << std::endl;
}


/**
 * @brief 写出pFind.cfg参数文件
 * 
 * @param cfg_path 写出路径
 */
void PFindParams::Write(const std::string &out_path) {
    FILE* cfg_file = OpenFileWrite(out_path);
    setvbuf ( cfg_file , NULL , _IOFBF , 102400 );
    std::string title = "# This is a standard pFind configure file\n# For help: mail to chihao@ict.ac.cn\n# Time: 2020-03-13 0:09:04\n\n[Version]\npFind_Version=EVA.3.0.11\n";
    std::string title_param = "\n[param]\n", title_filter = "\n[filter]\n", title_engine = "\n[engine]\n", title_file = "\n[file]\n", title_datalist = "\n[datalist]\n", title_quant = "\n[quant]\n", title_system = "\n[system]\n";
    fwrite(title.c_str(), title.length(), 1, cfg_file);
    fwrite(title_param.c_str(), title_param.length(), 1, cfg_file);
    auto write_line = [&] (auto one_param, const std::string &param_name) { // 写一行参数
		fwrite(param_name.c_str(), param_name.length(), 1, cfg_file);
        fwrite(std::to_string(one_param).c_str(), std::to_string(one_param).length(), 1, cfg_file);
        fwrite("\n", 1, 1, cfg_file);
    };
    auto write_line_string = [&] (const std::string & one_param, const std::string &param_name) { // 写一行参数
		fwrite(param_name.c_str(), param_name.length(), 1, cfg_file);
        fwrite(one_param.c_str(), one_param.length(), 1, cfg_file);
        fwrite("\n", 1, 1, cfg_file);
    };
    write_line(this->multip, "multip=");
    write_line(this->thread, "thread="); 
    write_line_string(this->activation_type, "activation_type="); 
    write_line(this->mstol, "mstol=");
    write_line(this->mstolppm, "mstolppm="); 
    write_line(this->msmstol, "msmstol="); 
    write_line(this->msmstolppm, "msmstolppm="); 
    write_line(this->temppepnum, "temppepnum="); 
    write_line(this->pepnum, "pepnum=");
    write_line(this->selectpeak, "selectpeak="); 
    write_line(this->maxprolen, "maxprolen="); 
    write_line(this->maxspec, "maxspec="); 
    write_line(this->IeqL, "IeqL="); 
    write_line(this->npep, "npep="); 
    write_line(this->maxdelta, "maxdelta="); 
    write_line_string(this->selectmod, "selectmod="); 
    write_line_string(this->fixmod, "fixmod="); 
    write_line(this->maxmod, "maxmod="); 
    write_line_string(this->enzyme, "enzyme="); 
    write_line(this->digest, "digest="); 
    write_line(this->max_clv_sites, "max_clv_sites=");
    fwrite(title_filter.c_str(), title_filter.length(), 1, cfg_file);
    write_line(this->psm_fdr, "psm_fdr=");
    write_line(this->psm_fdr_type, "psm_fdr_type="); 
    write_line(this->mass_lower, "mass_lower=");
    write_line(this->mass_upper, "mass_upper="); 
    write_line(this->len_lower, "len_lower=");
    write_line(this->len_upper, "len_upper="); 
    write_line(this->pep_per_pro, "pep_per_pro=");
    write_line(this->pro_fdr, "pro_fdr="); 
    fwrite(title_engine.c_str(), title_engine.length(), 1, cfg_file);
    write_line(this->open, "open=");
    write_line(this->open_tag_len, "open_tag_len="); 
    write_line(this->rest_tag_iteration, "rest_tag_iteration=");
    write_line(this->rest_tag_len, "rest_tag_len="); 
    write_line(this->rest_mod_num, "rest_mod_num=");
    write_line(this->salvo_iteration, "salvo_iteration="); 
    write_line(this->salvo_mod_num, "salvo_mod_num=");
    fwrite(title_file.c_str(), title_file.length(), 1, cfg_file);
    write_line_string(this->modpath, "modpath="); 
    write_line_string(this->fastapath, "fastapath=");
    write_line(this->auto_rev, "auto_rev="); 
    write_line_string(this->outputpath, "outputpath="); 
    fwrite(title_datalist.c_str(), title_datalist.length(), 1, cfg_file);
    write_line(this->msmsnum, "msmsnum="); 
    if (this->msmsfolder > 0)
        write_line(this->msmsfolder, "msmsfolder=");
    int msmscnt = 1;
    for (auto path = this->msmspaths.begin(), end_ = this->msmspaths.end(); path != end_; path++) {
        write_line_string(Strip(&(*path), "/\\"), "msmspath" + std::to_string(msmscnt) + "=");
        msmscnt++;
    }
    write_line_string(this->msmstype, "msmstype="); 
    fwrite(title_quant.c_str(), title_quant.length(), 1, cfg_file);
    write_line_string(this->quant, "quant="); 
    fwrite(title_system.c_str(), title_system.length(), 1, cfg_file);
    write_line_string(this->log, "log="); 
    fclose(cfg_file);
}

