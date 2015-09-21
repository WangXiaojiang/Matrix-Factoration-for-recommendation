#pragma once
#include<vector>
#include<map>
#include<string>
#include<Eigen/Eigen>
using namespace std;
using namespace Eigen;
const string occupations[21] = { "administrator",
"artist",
"doctor",
"educator",
"engineer",
"entertainment",
"executive",
"healthcare",
"homemaker",
"lawyer",
"librarian",
"marketing",
"none",
"other",
"programmer",
"retired",
"salesman",
"scientist",
"student",
"technician",
"writer"
};
const string sexs[2] = { "F", "M" };

class user
{
public:
	int id;
	int age;
	string sex;
	string occupation;
	VectorXd factor;
	double bias;
	double mean;
	double alpha;  //相似度权重系数
	map<int, int> rating;
	map<int, double> simS;
	map<int, double> simA;
	map<int, double> sim;
	map<int, double> topN;
public:
	user(int idx, int agex, string sexx, string occupationx, int factorNumber,int N)
	{
		id = idx;
		age = agex;
		sex = sexx;
		occupation = occupationx;
		factor = VectorXd::Random(factorNumber);
		bias = rand() % 5;
		mean = 0.0;
	}
	user(int idx, int agex, string sexx, string occupationx, int factorNumber)
	{
		id = idx;
		age = agex;
		sex = sexx;
		occupation = occupationx;
		factor = VectorXd::Random(factorNumber);
		bias = rand() % 5;
		mean = 0.0;
	}
	user(int idx, int agex, string sexx, string occupationx)
	{
		id = idx;
		age = agex;
		sex = sexx;
		occupation = occupationx;
	}
	user(){}
};
class item
{
public:
	int id;
	int catigories[19];
	VectorXd factor;
	double bias;
	double mean;
	double popular;
	double alpha;  //相似度权重系数
	map<int, int> rating;
	map<int, double> simS;
	map<int, double> simA;
	map<int, double> sim;
	map<int, double> topN;
public:
	item(int idx, int catigoriesx[19], int factorNumber,int N)
	{
		id = idx;
		for (int i = 0; i < 19; i++)
			catigories[i] = catigoriesx[i];
		factor = VectorXd::Random(factorNumber);
		bias = rand() % 5;
		mean = 0.0;
		popular = 0.0;
	}
	item(int idx, int catigoriesx[19], int factorNumber)
	{
		id = idx;
		for (int i = 0; i < 19; i++)
			catigories[i] = catigoriesx[i];
		factor = VectorXd::Random(factorNumber);
		bias = rand() % 5;
		mean = 0.0;
		popular = 0.0;
	}
	item(int idx, int catigoriesx[19])
	{
		id = idx;
		for (int i = 0; i < 19; i++)
			catigories[i] = catigoriesx[i];
	}
	item(){}
};

class age
{
public:
	int ages[4];
	map<int, VectorXd> factor;
	age(int factorNumber)
	{
		for (int i = 0; i < 4; i++)
		{
			ages[i] = i;
			VectorXd factorx = VectorXd::Random(factorNumber);
			factor.insert(pair<int, VectorXd>(ages[i], factorx));
		}
	}
	age()
	{}
};
class sex
{
public:
	int sexx[2];
	map<int, VectorXd> factor;
	map<string, int> nameToId;
	sex(int factorNumber)
	{
		for (int i = 0; i < 2; i++)
		{
			sexx[i] = i;
			VectorXd factorx = VectorXd::Random(factorNumber);
			nameToId.insert(pair<string, int>(sexs[i], sexx[i]));
			factor.insert(pair<int, VectorXd>(sexx[i], factorx));
		}
	}
	sex()
	{}
};
class occupation
{
public:
	int occupationx[21];
	map<int, VectorXd> factor;
	map<string, int> nameToId;
	occupation(int factorNumber)
	{
		for (int i = 0; i < 21; i++)
		{
			occupationx[i] = i;
			nameToId.insert(pair<string, int>(occupations[i], occupationx[i]));
			VectorXd factorx = VectorXd::Random(factorNumber);
			factor.insert(pair<int, VectorXd>(occupationx[i], factorx));
		}
	}
	occupation()
	{}

};
class catigries
{
public:
	int catigriesx[19];
	double factor[19];
	double popular[19];
	catigries()
	{
		for (int i = 0; i < 19; i++)
		{
			catigriesx[i] = i;
			factor[i] = double(rand() % 10)/10.0;
			popular[i] = 0.0;
		}
	}
};
/*
class matrixFactoration
{
public:
	string trainFileName;
	string testFileName;
	vector<vector<int>> trainRateData;
	vector<vector<int>> testRateData;
	map<int, user> userInfo;
	map<int, item> itemInfo;
	map<int, attribution> attributionInfo;
	matrixFactoration()
	{}
	int loadTrainFile(string fileName,string patten)
	{
		ifstream tf;
		string temp;
		tf.open(fileName);
		int flag = -1;
		while (getline(tf, temp))
		{
			flag = 1;
			vector<int> rate;
			vector<string> result = split(temp, patten);
			rate.push_back(std::atoi(result[0].c_str()));
			rate.push_back(std::atoi(result[1].c_str()));
			rate.push_back(std::atoi(result[2].c_str()));
			trainRateData.push_back(rate);
		}
		tf.close();
		return flag;
	}
	int loadTestFile(string fileName, string patten)
	{
		ifstream tf;
		string temp;
		tf.open(fileName);
		int flag = -1;
		while (getline(tf, temp))
		{
			flag = 1;
			vector<int> rate;
			vector<string> result = split(temp, patten);
			rate.push_back(std::atoi(result[0].c_str()));
			rate.push_back(std::atoi(result[1].c_str()));
			rate.push_back(std::atoi(result[2].c_str()));
			testRateData.push_back(rate);
		}
		tf.close();
		return flag;
	}
	int loadUserFile(string fileName, string patten,int attributionNumber=3)
	{
		ifstream tf;
		string temp;
		tf.open(fileName);
		int flag = -1;
		while (getline(tf, temp))
		{
			flag = 1;
			vector<string> attribution;
			vector<string> result = split(temp, patten);
			if (attributionNumber > result.size() - 1)
			{
				cout << "attribution number erro!" << endl;
				return -1;
			}			
			for (int i = 0; i < attributionNumber+1; i++)
			{
				attribution.push_back(result[i]);
			}
			userInfo.insert(pair<int,user>(atoi(result[0].c_str()),user(attribution)));
		}
		tf.close();
		if (!flag)
			cout << "Open "<<fileName<<" erro!"<< endl;
		return flag;
	}
};
*/
class stasticInfo
{
public:
	map<int, map<int, int>> age;
	map<int, map<int,double>> ageAvg;
	map<string, map<int, int>> sex;
	map<string, map<int, double>> sexAvg;
	map<string, map<int, int>> occupation;
	map<string, map<int, double>> occupationAvg;
	stasticInfo(map<int, int> itemInfo, map<int, double> itemInfoRate)
	{
		for (int i = 0; i < 5; i++)
		{
			age.insert(pair<int,map<int, int>>(i, itemInfo));
			ageAvg.insert(pair<int, map<int, double>>(i, itemInfoRate));
		}
		for (int i = 0; i < 2; i++)
		{
			sex.insert(pair<string, map<int, int>>(sexs[i], itemInfo));
			sexAvg.insert(pair<string, map<int, double>>(sexs[i], itemInfoRate));
		}
		for (int i = 0; i < 21; i++)
		{
			occupation.insert(pair<string, map<int, int>>(occupations[i], itemInfo));
			occupationAvg.insert(pair<string, map<int, double>>(occupations[i], itemInfoRate));
		}
	}
};