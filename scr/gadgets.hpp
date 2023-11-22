//
//  gadgets.hpp
//  gctb
//
//  Created by Jian Zeng on 14/06/2016.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#ifndef toolbox_hpp
#define toolbox_hpp

#include <sys/time.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <string>
#include <sys/stat.h>
#include <map>
#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;

namespace Gadget {

class Timer {
    time_t prev, curr;
public:
    Timer(){
        setTime();
    };
    void setTime(void);
    time_t getTime(void);
    time_t getElapse(void);
    string format(const time_t time);
    string getDate(void);
    void printElapse(void);
};

class Tokenizer : public vector<string> {
    // adopted from matvec
public:
    void getTokens(const string &str, const string &sep);
    int  getIndex(const string &str);
};

template <class T> class Recoder : public map<T,unsigned> {
    // adopted from matvec
    unsigned count;
public:
    Recoder(void){count=0;}
    unsigned code(T s){
        typename map<T,unsigned>::iterator mapit = this->find(s);
        if(mapit == this->end()){
            (*this)[s] = ++count;
            return count;
        }
        else {
            return (*mapit).second;
        }
    }
    void display_codes(ostream & os = cout){
        typename Recoder::iterator it;
        for (it=this->begin(); it!=this->end();it++){
            os << (*it).first << " " << (*it).second << endl;
        }
    }
};

// file process functions
string getFileName(const string &file);
string getFileSuffix(const string &file);
void fileExist(const string &filename);
bool directoryExist(const string &dirname);
bool createDirectory(const string& dirname);

// statistics functions
float calcMean(const VectorXf &vec);
float calcVariance(const VectorXf &vec);
float calcCovariance(const VectorXf &vec1, const VectorXf &vec2);
float calcCorrelation(const VectorXf &vec1, const VectorXf &vec2);
float calcRegression(const VectorXf &y, const VectorXf &x);
float findMedian(const VectorXf &vec);

}

#endif /* toolbox_hpp */
