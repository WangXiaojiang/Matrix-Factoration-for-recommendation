#include "approach.cpp"
#include "classObject.h"
#include "similarity.h"
const int FACTORNUMBER = 20;
double LEARNRATE = 0.01;
const double hp1 = 0.1;
const double hp2 = 0.1;
const double hp3 = 0.002;
const double hp4 = 0.1;
const int ITRMAX = 10000;
const int N = 10;
double AVG = 0.0;
const double thread = 0.7;
clock_t start, finish;
using namespace std;
void main()
{
	string userInfoFile = "u.user";
	string itemInfoFile = "u.item";
	string trainFile = "u3.base";
	string testFile = "u3.test";
	string userSimFile = "u3.sim";
	string itemSimFile = "i3.sim";
	map<int, user> userInfo;
	map<int, item> itemInfo;
	map<int, user>::iterator itUser;
	map<int, item>::iterator ititem;
	age age_case(FACTORNUMBER);
	sex sex_case(FACTORNUMBER);
	occupation occupation_case(FACTORNUMBER);
	ifstream userInfoFilePtr;
	userInfoFilePtr.open(userInfoFile);
	ifstream itemInfoFilePtr;
	itemInfoFilePtr.open(itemInfoFile);
	string temp;
	std::cout << "Start reading userInfoFile" << endl;
	for (; getline(userInfoFilePtr, temp);)
	{

		vector<string> temp1 = split(temp, "|");
		int idx = std::atoi(temp1[0].c_str());
		user user_case(std::atoi(temp1[0].c_str()), std::atoi(temp1[1].c_str()), temp1[2], temp1[3], FACTORNUMBER);
		userInfo.insert(pair<int, user>(idx, user_case));
	}
	userInfoFilePtr.close();
	std::cout << "UserInfoFile read over,start reading itemInfoFile" << endl;
	for (; getline(itemInfoFilePtr, temp);)
	{
		vector<string> temp1 = split(temp, "|");
		int idx = std::atoi(temp1[0].c_str());
		int catigariesx[19];
		for (int k = 0; k < 19; k++)
		{
			catigariesx[k] = std::atoi(temp1[4 + k].c_str());
		}
		item item_case(std::atoi(temp1[0].c_str()), catigariesx, FACTORNUMBER);
		itemInfo.insert(pair<int, item>(idx, item_case));
	}
	itemInfoFilePtr.close();
	std::cout << "ItemInfoFile read over,start read train set" << endl;

	ifstream trainFilePtr, testFilePtr;
	vector<vector<int>> train, test;
	trainFilePtr.open(trainFile);
	for (; getline(trainFilePtr, temp);)
	{
		vector<string> temp1 = split(temp, "\t");
		vector<int> ratex;
		int uidx = std::atoi(temp1[0].c_str());
		int iidx = std::atoi(temp1[1].c_str());
		int rate = std::atoi(temp1[2].c_str());
		ratex.push_back(uidx);
		ratex.push_back(iidx);
		ratex.push_back(rate);
		train.push_back(ratex);
		userInfo[uidx].rating[iidx] = rate;
		itemInfo[iidx].rating[uidx] = rate;
		AVG += rate;
	}
	trainFilePtr.close();
	AVG /= train.size();
	std::cout << "Train data read over,start read test data" << endl;
	testFilePtr.open(testFile);
	for (; getline(testFilePtr, temp);)
	{
		vector<string> temp1 = split(temp, "\t");
		vector<int> ratex;
		int uidx = std::atoi(temp1[0].c_str());
		int iidx = std::atoi(temp1[1].c_str());
		int rate = std::atoi(temp1[2].c_str());
		ratex.push_back(uidx);
		ratex.push_back(iidx);
		ratex.push_back(rate);
		test.push_back(ratex);
	}
	testFilePtr.close();
	std::cout << "Test data read over, star calculate user similarity" << endl;
	
	if (0)
	{
		ofstream userSimFilePtr, itemSimFilePtr;
		userSimFilePtr.open(userSimFile);
		itemSimFilePtr.open(itemSimFile);
		for (int i = 1; i <userInfo.size(); i++)
		{
			clock_t start, finish;
			start = clock();
			for (int j = i + 1; j <= userInfo.size(); j++)
			{
				userInfo[i].simS[j] = userInfo[j].simS[i] = similarity(userInfo[i].rating, userInfo[j].rating);
				double age = double(min(userInfo[i].age, userInfo[j].age)) / max(userInfo[i].age, userInfo[j].age);
				double gender=0,occupation=0;
				if (userInfo[i].sex == userInfo[j].sex)
					gender = 1;
				if (userInfo[i].occupation == userInfo[j].occupation)
					occupation = 1;			
				userInfo[i].simA[j] = userInfo[j].simA[i] =  (age + gender + occupation)/3;		
			}
			finish = clock();
			double dur = (double)(finish - start) / CLOCKS_PER_SEC;
			std::cout << "loop:" << i << "	total time : " << dur << endl;
		}
		std::cout << "User similarity calculate over,star calculate item similarity" << endl;
		for (int i = 1; i < itemInfo.size(); i++)
		{
			start = clock();
			for (int j = i + 1; j <= itemInfo.size(); j++)
			{
				itemInfo[i].simS[j] = itemInfo[j].simS[i] = similarity(itemInfo[i].rating, itemInfo[j].rating);
				itemInfo[i].simA[j] = itemInfo[j].simA[i] = vectorP(itemInfo[i].catigories, itemInfo[j].catigories, 19);
			}
			finish = clock();
			double dur = (double)(finish - start) / CLOCKS_PER_SEC;
			std::cout << "loop:" << i << "	total time : " << dur << endl;
		}
		std::cout << "Find top-N similar users" << endl;
		for (int i = 1; i <= userInfo.size(); i++)
		{
			start = clock();
			double alpha, beta;
			double number = userInfo[i].rating.size();
			if (number > 180)
			{
				 alpha = userInfo[i].alpha = 0.9;
				 beta = 0.1;
			}
			else
			{
				alpha = userInfo[i].alpha = number / 180.0;
				beta = 1 - alpha;
			}
			for (int j = 1; j <= userInfo.size(); j++)
			{
				if (i != j)
				{
					userInfo[i].sim[j] = alpha * userInfo[i].simS[j] + beta * userInfo[i].simA[j];
				}
			}			
			userSimFilePtr << i;
			map<int, double>::iterator itr = userInfo[i].sim.begin();
			while (itr != userInfo[i].sim.end())
			{
				userSimFilePtr <<"\t"<<itr->first<<":"<< itr->second;
				itr++;
			}
			userSimFilePtr <<endl;
			userInfo[i].topN = findTopN(userInfo[i], N, thread);
			finish = clock();
			double dur = (double)(finish - start) / CLOCKS_PER_SEC;
			std::cout << "Step : Find top-N similar users!" << "loop:" << i << "	total time : " << dur << endl;
		}
		std::cout << "Find top-N similar items" << endl;
		for (int i = 1; i <= itemInfo.size(); i++)
		{
			start = clock();
			double alpha, beta;
			double number =itemInfo[i].rating.size();
			if (number > 120)
			{
				alpha = itemInfo[i].alpha = 0.9;
				beta = 0.1;
			}
			else
			{
				alpha = itemInfo[i].alpha = number / 120.0;
				beta = 1 - alpha;
			}
			for (int j = 1; j <= itemInfo.size(); j++)
			{
				if (i != j)
				{
					itemInfo[i].sim[j] = alpha * itemInfo[i].simS[j] + beta * itemInfo[i].simA[j];
				}
			}			
			itemSimFilePtr << i;
			map<int, double>::iterator itr = itemInfo[i].sim.begin();
			while (itr != itemInfo[i].sim.end())
			{
				itemSimFilePtr <<"\t" << itr->first<< ":"<< itr->second;
				itr++;
			}
			itemSimFilePtr <<endl;
			itemInfo[i].topN = findTopN(itemInfo[i], N, thread);
			finish = clock();
			double dur = (double)(finish - start) / CLOCKS_PER_SEC;
			std::cout << "Step : Find top-N similar items!"<< "loop:" << i << "	total time : " << dur << endl;
		}
		userSimFilePtr.close();
		itemSimFilePtr.close();
	}
	else
	{
		std::cout << "Read similar file" << endl;
		map< int, map < int, double >> userSim = readSimFile(userSimFile);
		map<int, map<int, double>> itemSim = readSimFile(itemSimFile);
		for (int i = 1; i < userInfo.size(); i++)
		{
			userInfo[i].sim = userSim[i];
			userInfo[i].topN = findTopN(userInfo[i], N, thread);
		}
		for (int i = 1; i <itemInfo.size(); i++)
		{
			itemInfo[i].sim = itemSim[i];
			itemInfo[i].topN = findTopN(itemInfo[i], N, thread);
		}
	}
	std::cout << "Main function start!" << endl;
	//	GlobalMean(train, test);
	//	UserMean(userInfo, train, test);
	//	ItemMean(itemInfo, train, test);
	//	BiasMF(userInfo, itemInfo, train, test,AVG);
	//	BiasNAMF(userInfo, itemInfo, train, test, age_case, AVG);
	//	BiasNASOPMF(userInfo, itemInfo, train, test, age_case, sex_case, occupation_case, AVG);
	//	BiasNASMF(userInfo, itemInfo, train, test, age_case, sex_case, AVG,0.01,hp1,hp2,hp3,hp4,FACTORNUMBER,ITRMAX);
	//	SocialMF(userInfo,itemInfo,train,test,0.01,hp1,hp3,20,10000);
	//	RSVD(userInfo,itemInfo,train,test);
	BiasNSOMF(userInfo, itemInfo, train, test, age_case, sex_case,occupation_case, AVG, 0.01, hp1, hp2, hp3, hp4, FACTORNUMBER, ITRMAX);
}
