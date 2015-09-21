/*#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<map>
#include<cmath>
#include "classObject.h"
#include "similarity.h"
#include <time.h>
clock_t start, finish;
string userfile = "u.user";
string itemfile = "u.item";
string trainfile = "cs.train";
string testfile = "cs.test";
using namespace std;
void main()
{
	map<int, user> userInfo;
	//	map<int, item> itemInfo;
	map<int, user>::iterator itUser;
	//	map<int, item>::iterator ititem;
	ifstream userInfoFilePtr;
	userInfoFilePtr.open(userfile);
	//	ifstream itemInfoFilePtr;
	//	itemInfoFilePtr.open(itemfile);
	string temp;
	std::cout << "start reading userInfoFile" << endl;
	for (; getline(userInfoFilePtr, temp);)
	{
		vector<string> temp1 = split(temp, "|");
		int idx = std::atoi(temp1[0].c_str());
		user user_case(std::atoi(temp1[0].c_str()), std::atoi(temp1[1].c_str()), temp1[2], temp1[3]);
		userInfo.insert(pair<int, user>(idx, user_case));
	}
	userInfoFilePtr.close();
	ifstream trainFilePtr;
	trainFilePtr.open(trainfile);
	ifstream testFilePtr;
	testFilePtr.open(testfile);
	map<int, int> itemInfo;
	map<int, double> itemInfoRate;
	for (; getline(trainFilePtr, temp);)
	{
		vector<string> temp1 = split(temp, "\t");
		int uidx = std::atoi(temp1[0].c_str());
		int iidx = std::atoi(temp1[1].c_str());
		itemInfo.insert(pair<int, int>(iidx, 0));
		itemInfoRate.insert(pair<int, double>(iidx, 0.0));
	}
	trainFilePtr.close();
	for (; getline(testFilePtr, temp);)
	{
		vector<string> temp1 = split(temp, "\t");
		int uidx = std::atoi(temp1[0].c_str());
		int iidx = std::atoi(temp1[1].c_str());
		itemInfo.insert(pair<int, int>(iidx, 0));
		itemInfoRate.insert(pair<int, double>(iidx, 0.0));
	}
	testFilePtr.close();
	trainFilePtr.open(trainfile);
	stasticInfo stasticResult(itemInfo, itemInfoRate);
	for (; getline(trainFilePtr, temp);)
	{
		vector<string> temp1 = split(temp, "\t");
		int uidx = std::atoi(temp1[0].c_str());
		int iidx = std::atoi(temp1[1].c_str());
		int rate = std::atoi(temp1[2].c_str());
		stasticResult.age[userInfo[uidx].age / 20][iidx] += 1;
		stasticResult.ageAvg[userInfo[uidx].age / 20][iidx] = ((stasticResult.age[userInfo[uidx].age / 20][iidx] - 1)*stasticResult.ageAvg[userInfo[uidx].age / 20][iidx] + rate) / stasticResult.age[userInfo[uidx].age / 20][iidx];
		stasticResult.sex[userInfo[uidx].sex][iidx] += 1;
		stasticResult.sexAvg[userInfo[uidx].sex][iidx] = ((stasticResult.sex[userInfo[uidx].sex][iidx] - 1)*stasticResult.sexAvg[userInfo[uidx].sex][iidx] + rate) / stasticResult.sex[userInfo[uidx].sex][iidx];
		stasticResult.occupation[userInfo[uidx].occupation][iidx] += 1;
		stasticResult.occupationAvg[userInfo[uidx].occupation][iidx] = ((stasticResult.occupation[userInfo[uidx].occupation][iidx] - 1)*stasticResult.occupationAvg[userInfo[uidx].occupation][iidx] + rate) / stasticResult.occupation[userInfo[uidx].occupation][iidx];
	}
	trainFilePtr.close();
	testFilePtr.open(testfile);
	for (; getline(testFilePtr, temp);)
	{
		vector<string> temp1 = split(temp, "\t");
		int uidx = std::atoi(temp1[0].c_str());
		int iidx = std::atoi(temp1[1].c_str());
		int rate = std::atoi(temp1[2].c_str());
		stasticResult.age[userInfo[uidx].age / 20][iidx] += 1;
		stasticResult.ageAvg[userInfo[uidx].age / 20][iidx] = ((stasticResult.age[userInfo[uidx].age / 20][iidx] - 1)*stasticResult.ageAvg[userInfo[uidx].age / 20][iidx] + rate) / stasticResult.age[userInfo[uidx].age / 20][iidx];
		stasticResult.sex[userInfo[uidx].sex][iidx] += 1;
		stasticResult.sexAvg[userInfo[uidx].sex][iidx] = ((stasticResult.sex[userInfo[uidx].sex][iidx] - 1)*stasticResult.sexAvg[userInfo[uidx].sex][iidx] + rate) / stasticResult.sex[userInfo[uidx].sex][iidx];
		stasticResult.occupation[userInfo[uidx].occupation][iidx] += 1;
		stasticResult.occupationAvg[userInfo[uidx].occupation][iidx] = ((stasticResult.occupation[userInfo[uidx].occupation][iidx] - 1)*stasticResult.occupationAvg[userInfo[uidx].occupation][iidx] + rate) / stasticResult.occupation[userInfo[uidx].occupation][iidx];
	}
	testFilePtr.close();
	ofstream ageStastic[5], sexStastic[2], occupationStastic[21], ageStasticAvg[5], sexStasticAvg[2], occupationStasticAvg[21];
		for (int i = 0; i < 5; i++)
		{
			stringstream ss, sl;
			ss << "ageStastic" << i << ".txt";
			sl << "ageStasticAvg" << i << ".txt";
			string filename = ss.str();
			ageStastic[i].open(filename);
			ageStasticAvg[i].open(sl.str());
		}
		for (int i = 0; i < 5; i++)
		{
			map<int, int>::iterator itr = stasticResult.age[i].begin();
			map<int, double>::iterator itrf = stasticResult.ageAvg[i].begin();
			while (itr != stasticResult.age[i].end())
			{
				ageStastic[i] << itr->first << "\t" << itr->second << endl;
				itr++;
			}
			ageStastic[i].close();
			while (itrf != stasticResult.ageAvg[i].end())
			{
				ageStasticAvg[i] << itrf->first << "\t" << itrf->second << endl;
				itrf++;
			}
			ageStasticAvg[i].close();
		}
		for (int i = 0; i < 2; i++)
		{
			stringstream ss, sl;
			ss << "sexStastic" << i << ".txt";
			sl << "sexStasticAvg" << i << ".txt";
			string filename = ss.str();
			sexStastic[i].open(filename);
			sexStasticAvg[i].open(sl.str());
		}
		for (int i = 0; i < 2; i++)
		{
			map<int, int>::iterator itr = stasticResult.sex[sexs[i]].begin();
			map<int, double>::iterator itrf = stasticResult.sexAvg[sexs[i]].begin();
			while (itr != stasticResult.sex[sexs[i]].end())
			{
				sexStastic[i] << itr->first << "\t" << itr->second << endl;
				itr++;
			}
			sexStastic[i].close();
			while (itrf != stasticResult.sexAvg[sexs[i]].end())
			{
				sexStasticAvg[i] << itrf->first << "\t" << itrf->second << endl;
				itrf++;
			}
			sexStasticAvg[i].close();
		}
		for (int i = 0; i < 21; i++)
		{
			stringstream ss, sl;
			ss << "occupationStastic" << i << ".txt";
			sl << "occupationStasticAvg" << i << ".txt";
			string filename = ss.str();
			occupationStastic[i].open(filename);
			occupationStasticAvg[i].open(sl.str());
		}
		for (int i = 0; i < 21; i++)
		{
			map<int, int>::iterator itr = stasticResult.occupation[occupations[i]].begin();
			map<int, double>::iterator itrf = stasticResult.occupationAvg[occupations[i]].begin();
			while (itr != stasticResult.occupation[occupations[i]].end())
			{
				occupationStastic[i] << itr->first << "\t" << itr->second << endl;
				itr++;
			}
			occupationStastic[i].close();
			while (itrf != stasticResult.occupationAvg[occupations[i]].end())
			{
				occupationStasticAvg[i] << itrf->first << "\t" << itrf->second << endl;
				itrf++;
			}
			occupationStasticAvg[i].close();
		}
}
*/
