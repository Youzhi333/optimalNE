//  Algorithm CRM in optimal NE NeurIPS 2023 Cited: Youzhi Zhang, Bo An, V.S. Subrahmanian. Computing optimal Nash equilibria in multiplayer games. Proceedings of the 37th Conference on Neural Information Processing Systems (NeurIPS'23)
//  Change the parameters about numbers of players actions  
//  Games are randomly generated
//  Time limit is set as 1000s
//  Created by Youzhi Zhang .
//


 

#include <cassert>
#include "gurobi_c++.h"
#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <math.h>
#include <list>
#include <set>
#include <numeric>
#include <queue>
#include <stdexcept>
 
using namespace std;


 
int playerNumber = 5;
int actionNumber = 2;
vector<vector<int>> allCombinations;
float solvedTime = 0;
float CPUTime = 0;
void printcombos(const vector<vector<int>>& vec, vector<int>& index, int depth) {
    if (depth == index.size()) {
        vector<int> rA;
        for (int i = 0; i < depth; ++i) {
            rA.push_back(vec[i][index[i]]);
             
        }
        allCombinations.push_back(rA);
         
    }
    else {
        const vector<int>& myvec = vec[depth];
        long mylength = myvec.size();
        for (int i = 0; i < mylength; ++i) {
            index[depth] = i;
            printcombos(vec, index, depth + 1);
        }
    }
}
 

void binaryrecursivesetgeneration(vector<int> allTeamPlayers, vector<vector<int>>& allSet, map<int, vector<int>>& setlocation2childrenlocation, map<vector<int>, int>& set2location) {
 

    int level = ceil(log2(allTeamPlayers.size()));
    std::cout << "start!" << level << endl;
    if (pow(2, level) == allTeamPlayers.size()) {
        vector<vector<int>> lowerset;
        for (int p = 0; p < allTeamPlayers.size(); p++) {
            vector<int> smallset{ allTeamPlayers[p] };
            lowerset.push_back(smallset);
            allSet.push_back(smallset);
            set2location[smallset] = allSet.size() - 1;

        }

        vector<vector<int>> largeset(lowerset.size() / 2);
        for (int j = 0; j < lowerset.size() / 2; j++) {
            vector<int> left = lowerset[j * 2];
            vector<int> right = lowerset[j * 2 + 1];
            left.insert(left.end(), right.begin(), right.end());
            largeset[j] = left;
            allSet.push_back(left);
            set2location[left] = allSet.size() - 1;
             
        }
        lowerset = largeset;
        for (int i = 2; i <= level; i++) {
            vector<vector<int>> largeset(lowerset.size() / 2);
            for (int j = 0; j < lowerset.size() / 2; j++) {
                vector<int> left = lowerset[j * 2];
                vector<int> right = lowerset[j * 2 + 1];
                left.insert(left.end(), right.begin(), right.end());
                largeset[j] = left;
                allSet.push_back(left);
                set2location[left] = allSet.size() - 1;
                vector<int> lr{ set2location[lowerset[j * 2]] };
                lr.push_back(set2location[lowerset[j * 2 + 1]]);
                setlocation2childrenlocation[allSet.size() - 1] = lr;

            }
            lowerset = largeset;
        }

    }
    else {
        vector<int> set1;
        vector<int> set2;
        for (int p = 0; p < pow(2, level - 1); p++) {
            set1.push_back(allTeamPlayers[p]);
        }
        for (int p = pow(2, level - 1); p < allTeamPlayers.size(); p++) {
            set2.push_back(allTeamPlayers[p]);
        }

        if (3 * pow(2, level - 2) <= allTeamPlayers.size()) {


            vector<vector<int>> lowerset;
            for (int p = 0; p < set1.size(); p++) {
                vector<int> smallset{ set1[p] };
                lowerset.push_back(smallset);
                allSet.push_back(smallset);
                set2location[smallset] = allSet.size() - 1;

            }

            vector<vector<int>> largeset(lowerset.size() / 2);
            for (int j = 0; j < lowerset.size() / 2; j++) {
                vector<int> left = lowerset[j * 2];
                vector<int> right = lowerset[j * 2 + 1];
                left.insert(left.end(), right.begin(), right.end());
                largeset[j] = left;
                allSet.push_back(left);
                set2location[left] = allSet.size() - 1;
                
            }
            lowerset = largeset;
            for (int i = 2; i <= level - 1; i++) {//set1 has level-1
                vector<vector<int>> largeset(lowerset.size() / 2);
                for (int j = 0; j < lowerset.size() / 2; j++) {
                    vector<int> left = lowerset[j * 2];
                    vector<int> right = lowerset[j * 2 + 1];
                    left.insert(left.end(), right.begin(), right.end());
                    largeset[j] = left;
                    allSet.push_back(left);
                    set2location[left] = allSet.size() - 1;
                    vector<int> lr{ set2location[lowerset[j * 2]] };
                    lr.push_back(set2location[lowerset[j * 2 + 1]]);
                    setlocation2childrenlocation[allSet.size() - 1] = lr;

                }
                lowerset = largeset;
            }
            binaryrecursivesetgeneration(set2, allSet, setlocation2childrenlocation, set2location);
            allSet.push_back(allTeamPlayers);
            set2location[allTeamPlayers] = allSet.size() - 1;
            vector<int> lr{ set2location[set1] };
            lr.push_back(set2location[set2]);
            setlocation2childrenlocation[allSet.size() - 1] = lr;

        }
        else {
            vector<vector<int>> lowerset;
            for (int p = 0; p < set1.size(); p++) {
                vector<int> smallset{ set1[p] };
                lowerset.push_back(smallset);
                allSet.push_back(smallset);
                set2location[smallset] = allSet.size() - 1;

            }

            vector<vector<int>> largeset(lowerset.size() / 2);
            for (int j = 0; j < lowerset.size() / 2; j++) {
                vector<int> left = lowerset[j * 2];
                vector<int> right = lowerset[j * 2 + 1];
                left.insert(left.end(), right.begin(), right.end());
                largeset[j] = left;
                allSet.push_back(left);
                set2location[left] = allSet.size() - 1;
                
            }
            lowerset = largeset;
            for (int i = 2; i <= level - 2; i++) { 
                vector<vector<int>> largeset(lowerset.size() / 2);
                for (int j = 0; j < lowerset.size() / 2; j++) {
                    vector<int> left = lowerset[j * 2];
                    vector<int> right = lowerset[j * 2 + 1];
                    left.insert(left.end(), right.begin(), right.end());
                    largeset[j] = left;
                    allSet.push_back(left);
                    set2location[left] = allSet.size() - 1;
                    vector<int> lr{ set2location[lowerset[j * 2]] };
                    lr.push_back(set2location[lowerset[j * 2 + 1]]);
                    setlocation2childrenlocation[allSet.size() - 1] = lr;

                }
                lowerset = largeset;
            }
            binaryrecursivesetgeneration(set2, allSet, setlocation2childrenlocation, set2location);
            vector<int> left = lowerset[1];
            vector<int> right = set2;
            left.insert(left.end(), right.begin(), right.end());

            allSet.push_back(left);
            set2location[left] = allSet.size() - 1;
            vector<int> lr{ set2location[lowerset[1]] };
            lr.push_back(set2location[set2]);
            setlocation2childrenlocation[allSet.size() - 1] = lr;

           
            allSet.push_back(allTeamPlayers);
            set2location[allTeamPlayers] = allSet.size() - 1;
            vector<int> alllr{ set2location[lowerset[0]] };
            alllr.push_back(set2location[left]);
            setlocation2childrenlocation[allSet.size() - 1] = alllr;

        }

    }




}


 
void CRM() {

    //#define PLAYERNUMBER 3
    int PLAYERNUMBER = playerNumber;
    int player = PLAYERNUMBER - 1;



     
    float minU = -INFINITY;
    float maxU = INFINITY;
    vector<vector<vector<float>>> termToPayoffs;
    vector<int> numberOfActions(PLAYERNUMBER);

    vector<vector<int>> terms;

    vector<int> allTeamPlayers;
    for (int p = 0; p < PLAYERNUMBER; p++) {
        numberOfActions[p] = actionNumber;
        if (p != player) {
            allTeamPlayers.push_back(p);
        }
    }
    numberOfActions[player] = actionNumber;
    termToPayoffs.resize(numberOfActions[player]);



    if (PLAYERNUMBER == 3) {
        for (int i = 0; i < numberOfActions[allTeamPlayers[0]]; i++) {
            for (int j = 0; j < numberOfActions[allTeamPlayers[1]]; j++) {
                terms.push_back({ i,j });
            }
        }
    }
    if (PLAYERNUMBER == 4) {
        for (int i = 0; i < numberOfActions[allTeamPlayers[0]]; i++) {
            for (int j = 0; j < numberOfActions[allTeamPlayers[1]]; j++) {
                for (int k = 0; k < numberOfActions[allTeamPlayers[2]]; k++) {

                    terms.push_back({ i,j,k });
                }
            }
        }
    }
    if (PLAYERNUMBER == 5) {
        for (int i = 0; i < numberOfActions[allTeamPlayers[0]]; i++) {
            for (int j = 0; j < numberOfActions[allTeamPlayers[1]]; j++) {
                for (int k = 0; k < numberOfActions[allTeamPlayers[2]]; k++) {
                    for (int l = 0; l < numberOfActions[allTeamPlayers[3]]; l++) {
                        terms.push_back({ i,j,k,l });
                    }
                }
            }
        }
    }

    if (PLAYERNUMBER == 6) {
        for (int i = 0; i < numberOfActions[allTeamPlayers[0]]; i++) {
            for (int j = 0; j < numberOfActions[allTeamPlayers[1]]; j++) {
                for (int k = 0; k < numberOfActions[allTeamPlayers[2]]; k++) {
                    for (int l = 0; l < numberOfActions[allTeamPlayers[3]]; l++) {
                        for (int m = 0; m < numberOfActions[allTeamPlayers[4]]; m++) {
                            terms.push_back({ i,j,k,l,m });
                        }
                    }
                }
            }
        }
    }

    if (PLAYERNUMBER == 7) {
        for (int i = 0; i < numberOfActions[allTeamPlayers[0]]; i++) {
            for (int j = 0; j < numberOfActions[allTeamPlayers[1]]; j++) {
                for (int k = 0; k < numberOfActions[allTeamPlayers[2]]; k++) {
                    for (int l = 0; l < numberOfActions[allTeamPlayers[3]]; l++) {
                        for (int m = 0; m < numberOfActions[allTeamPlayers[4]]; m++) {
                            for (int n = 0; n < numberOfActions[allTeamPlayers[5]]; n++) {
                                terms.push_back({ i,j,k,l,m,n });
                            }
                        }
                    }
                }
            }
        }
    }

    if (PLAYERNUMBER == 8) {
        for (int i = 0; i < numberOfActions[allTeamPlayers[0]]; i++) {
            for (int j = 0; j < numberOfActions[allTeamPlayers[1]]; j++) {
                for (int k = 0; k < numberOfActions[allTeamPlayers[2]]; k++) {
                    for (int l = 0; l < numberOfActions[allTeamPlayers[3]]; l++) {
                        for (int m = 0; m < numberOfActions[allTeamPlayers[4]]; m++) {
                            for (int n = 0; n < numberOfActions[allTeamPlayers[5]]; n++) {
                                for (int a = 0; a < numberOfActions[allTeamPlayers[6]]; a++) {
                                    terms.push_back({ i,j,k,l,m,n,a });
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    if (PLAYERNUMBER == 9) {
        for (int i = 0; i < numberOfActions[allTeamPlayers[0]]; i++) {
            for (int j = 0; j < numberOfActions[allTeamPlayers[1]]; j++) {
                for (int k = 0; k < numberOfActions[allTeamPlayers[2]]; k++) {
                    for (int l = 0; l < numberOfActions[allTeamPlayers[3]]; l++) {
                        for (int m = 0; m < numberOfActions[allTeamPlayers[4]]; m++) {
                            for (int n = 0; n < numberOfActions[allTeamPlayers[5]]; n++) {
                                for (int a = 0; a < numberOfActions[allTeamPlayers[6]]; a++) {
                                    for (int b = 0; b < numberOfActions[allTeamPlayers[7]]; b++) {
                                        terms.push_back({ i,j,k,l,m,n,a,b });
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    if (PLAYERNUMBER == 10) {
        for (int i = 0; i < numberOfActions[allTeamPlayers[0]]; i++) {
            for (int j = 0; j < numberOfActions[allTeamPlayers[1]]; j++) {
                for (int k = 0; k < numberOfActions[allTeamPlayers[2]]; k++) {
                    for (int l = 0; l < numberOfActions[allTeamPlayers[3]]; l++) {
                        for (int m = 0; m < numberOfActions[allTeamPlayers[4]]; m++) {
                            for (int n = 0; n < numberOfActions[allTeamPlayers[5]]; n++) {
                                for (int a = 0; a < numberOfActions[allTeamPlayers[6]]; a++) {
                                    for (int b = 0; b < numberOfActions[allTeamPlayers[7]]; b++) {
                                        for (int c = 0; c < numberOfActions[allTeamPlayers[8]]; c++) {
                                            terms.push_back({ i,j,k,l,m,n,a,b,c });
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if (PLAYERNUMBER == 12) {
        for (int i = 0; i < numberOfActions[allTeamPlayers[0]]; i++) {
            for (int j = 0; j < numberOfActions[allTeamPlayers[1]]; j++) {
                for (int k = 0; k < numberOfActions[allTeamPlayers[2]]; k++) {
                    for (int l = 0; l < numberOfActions[allTeamPlayers[3]]; l++) {
                        for (int m = 0; m < numberOfActions[allTeamPlayers[4]]; m++) {
                            for (int n = 0; n < numberOfActions[allTeamPlayers[5]]; n++) {
                                for (int a = 0; a < numberOfActions[allTeamPlayers[6]]; a++) {
                                    for (int b = 0; b < numberOfActions[allTeamPlayers[7]]; b++) {
                                        for (int c = 0; c < numberOfActions[allTeamPlayers[8]]; c++) {
                                            for (int d = 0; d < numberOfActions[allTeamPlayers[9]]; d++) {
                                                for (int e = 0; e < numberOfActions[allTeamPlayers[10]]; e++) {
                                                    terms.push_back({ i,j,k,l,m,n,a,b,c,d,e });
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    
    for (int i = 0; i < numberOfActions[player]; i++) {
        for (int j = 0; j < terms.size(); j++) {
            vector<float> payoff(PLAYERNUMBER);
            for (int p = 0; p < PLAYERNUMBER; p++) {
                payoff[p] = std::rand() / ((RAND_MAX + 1u) / 100);
                
            }
             
            termToPayoffs[i].push_back(payoff);
             
        }
        

    }
     
    


     
 

 
    map<vector<int>, int> set2location;

    vector<vector<int>> allSet;
    map<int, vector<int>> setlocation2childrenlocation; 

    binaryrecursivesetgeneration(allTeamPlayers, allSet, setlocation2childrenlocation, set2location);

    for (int p = 0; p < allTeamPlayers.size(); p++) {


        vector<int> temp = allTeamPlayers;
        temp[p] = player;
        allSet.push_back(temp);
        set2location[temp] = allSet.size() - 1;
        vector<int> currentset = allTeamPlayers;
        bool shouldcontinue = true;
        if (currentset.size() < 3) {
            shouldcontinue = false;
        }
        int currentindex = allSet.size() - 1;
        while (shouldcontinue) {

            vector<int> children = setlocation2childrenlocation[set2location[currentset]];
            vector<int> leftchild = allSet[children[0]];
            vector<int> rightchild = allSet[children[1]];

            if (find(leftchild.begin(), leftchild.end(), allTeamPlayers[p]) != leftchild.end()) {
                currentset = leftchild;
                int index = distance(leftchild.begin(), find(leftchild.begin(), leftchild.end(), allTeamPlayers[p]));
                leftchild[index] = player;
                allSet.push_back(leftchild);
                set2location[leftchild] = allSet.size() - 1;
                vector<int> lr{ set2location[leftchild] };
                lr.push_back(set2location[rightchild]);
                setlocation2childrenlocation[currentindex] = lr;
                currentindex = allSet.size() - 1;


            }
            else {
                currentset = rightchild;
                int index = distance(rightchild.begin(), find(rightchild.begin(), rightchild.end(), allTeamPlayers[p]));
                rightchild[index] = player;
                allSet.push_back(rightchild);
                set2location[rightchild] = allSet.size() - 1;
                vector<int> lr{ set2location[leftchild] };
                lr.push_back(set2location[rightchild]);
                setlocation2childrenlocation[currentindex] = lr;
                currentindex = allSet.size() - 1;
            }


            if (currentset.size() < 3) {
                shouldcontinue = false;
            }
        }
    }
 


    vector<vector<vector<int>>> AllTerms(allSet.size());
    vector<map<vector<int>, int>> fromActionSetToIndex(allSet.size());
    for (int i = 0; i < allSet.size(); i++) {
        if (allSet[i].size() > 1) {
             

            vector<vector<int>> eachResourceAction;
            
            for (int p = 0; p < allSet[i].size(); p++) {
               
                vector<int> actionset;
                for (int j = 0; j < numberOfActions[allSet[i][p]]; j++) {
                    actionset.push_back(j);
                }
                eachResourceAction.push_back(actionset);
            }
           
            vector<int> actCombination;
            actCombination.resize(eachResourceAction.size());
            allCombinations.clear();
            printcombos(eachResourceAction, actCombination, 0);
            AllTerms[i] = allCombinations;
            for (int j = 0; j < allCombinations.size(); j++) {
                fromActionSetToIndex[i][allCombinations[j]] = j;
                 
            }
        }
        

    }

    /*****************************billinear term***************************************************************************/
    GRBEnv env = GRBEnv();

    GRBModel modelEQ = GRBModel(env);
    GRBVar** rSQ = new GRBVar * [PLAYERNUMBER];
    GRBVar** U = new GRBVar * [PLAYERNUMBER];
    GRBVar** B = new GRBVar * [PLAYERNUMBER];
 
    GRBVar** Vr = new GRBVar * [PLAYERNUMBER];
    GRBVar** Ur = new GRBVar * [PLAYERNUMBER];
    for (int p = 0; p < PLAYERNUMBER; p++) {
        rSQ[p] = new GRBVar[numberOfActions[p]];
        U[p] = new GRBVar[numberOfActions[p]];
        B[p] = new GRBVar[numberOfActions[p]];
        Ur[p] = new GRBVar[numberOfActions[p]];
        Vr[p] = new GRBVar[numberOfActions[p]];
        for (int i = 0; i < numberOfActions[p]; i++) {
            string aw = "rsq_" + std::to_string(p) + "_" + std::to_string(i);
            rSQ[p][i] = modelEQ.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, aw);

            U[p][i] = modelEQ.addVar(0, 100, 0.0, GRB_CONTINUOUS, "U_" + std::to_string(p) + "_" + std::to_string(i));
            Ur[p][i] = modelEQ.addVar(0, 100, 0.0, GRB_CONTINUOUS);
            Vr[p][i] = modelEQ.addVar(0, 100, 0.0, GRB_CONTINUOUS);
            B[p][i] = modelEQ.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        }
        GRBLinExpr linexpr;
         
        for (int j = 0; j < numberOfActions[p]; j++) {
            linexpr += rSQ[p][j];
        }
        modelEQ.addConstr(linexpr == 1);
        // }
    }
    

    GRBVar** w = new GRBVar * [AllTerms.size()];

  
    for (int i = 0; i < AllTerms.size(); i++) {

        w[i] = new GRBVar[AllTerms[i].size()];
        for (int j = 0; j < AllTerms[i].size(); j++) {
            string aw = "w_" + std::to_string(i) + "_" + std::to_string(j);
            w[i][j] = modelEQ.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, aw);
        }
    }
    for (int i = 0; i < AllTerms.size(); i++) {
        vector<int> playerVector = allSet[i];
        if (playerVector.size() == 2) {
            for (int bl = 0; bl < AllTerms[i].size(); bl++) {
                int ax = bl;
                string aw = "w_" + std::to_string(i) + "_" + std::to_string(ax);
                vector<int> term = AllTerms[i][bl];
                modelEQ.addQConstr(w[i][bl] == rSQ[playerVector[0]][term[0]] * rSQ[playerVector[1]][term[1]], aw);
            }
            for (int p = 0; p < playerVector.size(); p++) {
                for (int ai = 0; ai < numberOfActions[playerVector[p]]; ai++) {
                    GRBLinExpr constraints;
                    for (int bl = 0; bl < AllTerms[i].size(); bl++) {
                        vector<int> term = AllTerms[i][bl];
                        
                        if (ai == term[p]) {
                            
                            constraints += w[i][bl];
                        }
                    }
                    string aw = "wrsq_" + std::to_string(i) + "_" + std::to_string(playerVector[p]) + "_" + std::to_string(ai);
                    modelEQ.addConstr(constraints == rSQ[playerVector[p]][ai], aw);
                }
            }
        }
        else {
            if (playerVector.size() > 2) {
                vector<int> children = setlocation2childrenlocation[i];
                vector<int> leftchild = allSet[children[0]];
                vector<int> rightchild = allSet[children[1]];
                if (rightchild.size() == 1) {
                    for (int bl = 0; bl < AllTerms[i].size(); bl++) {
                        int ax = bl;
                        string aw = "w_" + std::to_string(i) + "_" + std::to_string(ax);
                        vector<int> term = AllTerms[i][bl];
                        vector<int> subterm;
                         
                        for (int p = 0; p < playerVector.size() - 1; p++) {
                            subterm.push_back(term[p]);
                             
                        }
                         
                        int index2left = fromActionSetToIndex[children[0]][subterm];
                        modelEQ.addQConstr(w[i][bl] == w[children[0]][index2left] * rSQ[playerVector[playerVector.size() - 1]][term[playerVector.size() - 1]], aw);
                    }

                }
                else {
                    for (int bl = 0; bl < AllTerms[i].size(); bl++) {
                        int ax = bl;
                        string aw = "w_" + std::to_string(i) + "_" + std::to_string(ax);
                        vector<int> term = AllTerms[i][bl];
                        vector<int> subterm;
                        vector<int> subterm2;
                        for (int p = 0; p < leftchild.size(); p++) {
                            subterm.push_back(term[p]);
                             
                        }
                        for (int p = leftchild.size(); p < playerVector.size(); p++) {
                            subterm2.push_back(term[p]);
                            
                        }
                        
                        int index2left = fromActionSetToIndex[children[0]][subterm];
                        int index2right = fromActionSetToIndex[children[1]][subterm2];
                        modelEQ.addQConstr(w[i][bl] == w[children[0]][index2left] * w[children[1]][index2right], aw);
                    }
                }
                 
                for (int p = 0; p < playerVector.size(); p++) {
                    for (int ai = 0; ai < numberOfActions[playerVector[p]]; ai++) {
                        GRBLinExpr constraints;
                        for (int bl = 0; bl < AllTerms[i].size(); bl++) {
                            vector<int> term = AllTerms[i][bl];
                            
                            if (ai == term[p]) {
                                 
                                constraints += w[i][bl];
                            }
                        }
                        string aw = "wrsq_" + std::to_string(i) + "_" + std::to_string(playerVector[p]) + "_" + std::to_string(ai);
                        modelEQ.addConstr(constraints == rSQ[playerVector[p]][ai], aw);
                    }
                }
            }


        }
        if (playerVector.size() >= 2) {
            GRBLinExpr constraints;
            for (int bl = 0; bl < AllTerms[i].size(); bl++) {
                constraints += w[i][bl];
            }
            modelEQ.addConstr(constraints == 1);
        }

    }


    for (int i = 0; i < AllTerms.size(); i++) {
        vector<int> playerVector = allSet[i];
        if (playerVector.size() > 2) {
            vector<int> children = setlocation2childrenlocation[i];
            vector<int> leftchild = allSet[children[0]];
            vector<int> rightchild = allSet[children[1]];
            if (rightchild.size() == 1) {
                for (int bl = 0; bl < AllTerms[children[0]].size(); bl++) {
                    GRBLinExpr constraints;
                    for (int j = 0; j < numberOfActions[playerVector[playerVector.size() - 1]]; j++) {
                        vector<int> subterm = AllTerms[children[0]][bl];
                        subterm.push_back(j);

                        int index2 = fromActionSetToIndex[i][subterm];
                        constraints += w[i][index2];
                    }
                    modelEQ.addConstr(constraints == w[children[0]][bl]);//the single elelment right part has been considered

                }


                


            }
            else {

                for (int bl = 0; bl < AllTerms[children[0]].size(); bl++) {
                    GRBLinExpr constraints;
                    for (int j = 0; j < AllTerms[children[1]].size(); j++) {
                        vector<int> subterm = AllTerms[children[0]][bl];
                        subterm.insert(subterm.end(), AllTerms[children[1]][j].begin(), AllTerms[children[1]][j].end());

                        int index2 = fromActionSetToIndex[i][subterm];
                        constraints += w[i][index2];
                    }
                    modelEQ.addConstr(constraints == w[children[0]][bl]);
                }

                for (int bl = 0; bl < AllTerms[children[1]].size(); bl++) {
                    GRBLinExpr constraints;
                    for (int j = 0; j < AllTerms[children[0]].size(); j++) {
                        vector<int> subterm = AllTerms[children[0]][j];
                        subterm.insert(subterm.end(), AllTerms[children[1]][bl].begin(), AllTerms[children[1]][bl].end());

                        int index2 = fromActionSetToIndex[i][subterm];
                        constraints += w[i][index2];
                    }
                    modelEQ.addConstr(constraints == w[children[1]][bl]);
                }



            }

        }

    }

 



    GRBLinExpr obj;
    GRBVar* Vobj = new GRBVar[PLAYERNUMBER];

    for (int p = 0; p < PLAYERNUMBER; p++) {
        Vobj[p] = modelEQ.addVar(0, 100, 0.0, GRB_CONTINUOUS);
        if (p == player) {
            obj += Vobj[p];
            
        }
    }
 
    modelEQ.setObjective(obj, GRB_MAXIMIZE);



    GRBQuadExpr constraints;
    GRBQuadExpr constraints3;

    for (int i = 0; i < numberOfActions[player]; i++) {
        GRBLinExpr QExpr;
        for (int j = 0; j < terms.size(); j++) {
 
            int index1 = set2location[allTeamPlayers];
 
            QExpr += termToPayoffs[i][j][player] * w[index1][j]; 
        }

        modelEQ.addConstr(QExpr == U[player][i]);
 
        modelEQ.addConstr(rSQ[player][i] <= 1 - B[player][i]);
        modelEQ.addConstr(Vobj[player] - U[player][i] <= 1000 * B[player][i]);

 
        modelEQ.addConstr(QExpr <= Vobj[player]);
    }
 

    for (int p = 0; p < allTeamPlayers.size(); p++) {
        GRBQuadExpr constraints2;
        GRBQuadExpr constraints3;
        vector<int> playerSet = allTeamPlayers;
        playerSet[p] = player;
 
        int index1 = set2location[playerSet];

        for (int i = 0; i < numberOfActions[allTeamPlayers[p]]; i++) {
            GRBLinExpr constraints;


            for (int bl = 0; bl < AllTerms[index1].size(); bl++) {
                vector<int> term = AllTerms[index1][bl];
                int playeraction = term[p];
                vector<int> newTerm = term;
                newTerm[p] = i;
 
                int index2 = set2location[allTeamPlayers];
                int index3 = fromActionSetToIndex[index2][newTerm];
 
                constraints += termToPayoffs[playeraction][index3][allTeamPlayers[p]] * w[index1][bl];
            }
            string aw = "w_" + std::to_string(p) + "_" + std::to_string(i) + "_";
            modelEQ.addConstr(constraints <= Vobj[allTeamPlayers[p]], aw + "V");
            modelEQ.addConstr(constraints == U[allTeamPlayers[p]][i], aw + "U");
 
            modelEQ.addConstr(rSQ[allTeamPlayers[p]][i] <= 1 - B[allTeamPlayers[p]][i]);
            modelEQ.addConstr(Vobj[allTeamPlayers[p]] - U[allTeamPlayers[p]][i] <= 1000 * B[allTeamPlayers[p]][i]);
 
        }
         
    }
 


    try {
 
        
        modelEQ.set(GRB_IntParam_NonConvex, 2);
       // modelEQ.set(GRB_DoubleParam_TimeLimit, 1000);
        modelEQ.optimize();

        

        
    }
    catch (GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
    catch (...) {
        cout << "Exception during optimization" << endl;
    }
    delete[] rSQ;



}

 

void RandomGame() {

    
     
    for (int index = 1; index <= 30; index++) {
        std::srand(20201125 + index * 10);
        
        CRM(); 
        

    }
     
     
}


 

int
main(int   argc,
    char* argv[])
{
     
    RandomGame();
    
    return 0;
}
