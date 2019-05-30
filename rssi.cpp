#include <iostream>
#include <armadillo>
#include <cmath>
#include <map>

using namespace arma;
using namespace std;

 
struct APList
{
    char   uuid;
    double x;        
    double y;
    double height;
    double rssi;     // receive one of aps rssi;
    double distance; // distance to ap (for test A ,n value);

    double A;
    double N;
    double weight;
};

double AP_A = 60.75;        // 1 meter Rssi Value;
double AP_n = 1.3026;        // interfence factors;
double totalWeight = 0;
map<int ,APList> APMap;
map<int ,APList> APCal;
mat result(2,1);

/*
 * d = 10^(abs(Rssi)- A)/()10*n)
 * ignore hight
 */
double GetDisFromRssi(APList aplist)
{
    double rawDis = 0.0;
    double power = (aplist.rssi - AP_A)/(10 * AP_n);
    rawDis = pow(10, power);

    return round(rawDis * 100)/100;
}

mat CaculateByAPList(APList apArray[])
{
    
    mat A(2,2);
    mat b(2,1);
    //matrix A init
    for(int i = 0; i < 2; i++)
    {
        A(i,0) = 2 * (apArray[i].x - apArray[2].x);
        A(i,1) = 2 * (apArray[i].y - apArray[2].y);
    }


    //matrix b init
    for(int i = 0; i < 2; i++)
    {
        double dTmp = 0.0;
        dTmp = pow(apArray[i].x, 2)
            - pow(apArray[2].x, 2)
            + pow(apArray[i].y, 2)
            - pow(apArray[2].y, 2)
            + pow(apArray[2].distance, 2)
            - pow(apArray[i].distance, 2);

        b(i,0) = dTmp;
    }
    //A.print("\nA=:");
    //b.print("b=:");
    //X=(A^T * A)^-1 * A^T * b

    double weight = 0;
    printf("AP:");
    for(int i = 0; i < 3; i++)
    {
        printf(" %d",apArray[i].uuid);
        weight += apArray[i].weight;
    }

    mat x(2,1);
    try {
        x = solve(A, b);
    }catch(runtime_error err){
        cout<< err.what() <<endl;
        x << 0 <<endr << 0 <<endr;
        weight = 0;
    }



    totalWeight += weight;

    x.print("; x =:");
    return x * weight;
}



//read aps info from file and set value;
int SetAPList(char *APFile, int flag = 0)
{
    FILE *fp = fopen(APFile, "r");
    if(fp == NULL)
    {
        printf("read test data file:[%s] failure!", APFile);
        return -1;
    }

    char line[1024];
    memset(line, 0x00 ,sizeof(line));
    APList ap;
    int rssi100;
    while(!feof(fp))
    {
        memset(&ap, 0x00, sizeof(ap));
        rssi100 = 0;
        fgets(line, 1024, fp);

        if(line[0] == '#')
            continue;

        //ap name or id
        char *pTmp = strtok(line, ",");
        if(!pTmp)
            continue;
        ap.uuid = atoi(pTmp);
        
        //ap location x
        pTmp = strtok(NULL, ",");
        if(!pTmp)
            continue;
        ap.x = atof(pTmp);
        //ap location y
        pTmp = strtok(NULL, ",");
        if(!pTmp)
            continue;
        ap.y = atof(pTmp);
        //ap  rssi
        pTmp = strtok(NULL, ",");
        if(!pTmp)
            continue;
        ap.rssi = abs(atof(pTmp));
        rssi100 =  100*ap.rssi;
        //ap  distance
        pTmp = strtok(NULL, ",");
        if(!pTmp)
            continue;

        if(!flag)
        {
            ap.distance = GetDisFromRssi(ap);
            printf("ap[%d]:x[%.2f] y[%.2f] rssi=%.2f, dis=%.2f\n",ap.uuid ,ap.x,ap.y,ap.rssi, ap.distance);
            if(ap.distance > 20)
                continue;
            ap.weight = 1.0/ap.distance;
            APMap.insert(make_pair(rssi100, ap));
        }
        else
        {
            ap.distance = atof(pTmp);
            APCal.insert(make_pair(rssi100, ap));
            if(ap.distance == 1)
            {
                AP_A = ap.rssi;
            }
        }
    }


    fclose(fp);

    return 0;
}

int CaculateAnByAPList(char* apFile)
{
    if(SetAPList(apFile, 1))
    {
        printf("[%s-%d].ERROR: calu fail for :Set APlist faile!\n", __FILE__,__LINE__);
        return 0;
    }

    double nTmp = 0;
    int n = 0;
    map<int, APList >::iterator it = APCal.begin();
    for(int i=0;i < APCal.size();i++,it++)
    {
        if(it->second.rssi != AP_A)
        {
            n++;
            it->second.N = (it->second.rssi - AP_A)/(10 * log10(it->second.distance));
        }
        nTmp += it->second.N; 

        printf("ap[%d],n[%f]\n",it->second.uuid, it->second.N);
    }
    AP_n = nTmp /n;


    return 0;
}
int combineNum(int n, int m)
{
    if(n == 0||m ==0)
        return 0;

    int a =1;
    int b =1;
    for(int i = 0; i < m; i++, n--)
        a = a * n;

    for(int i = 1; i < m+1; i++)
        b = b * i;

    return a/b;
}

int h = 1;
void combine(int n,int m,int *b,const int M)
{
    for(int j=n;j>=m;j--)
    {
        b[m-1]=j-1;
        if(m>1)
            combine(j-1,m-1,b,M);
        else
        {
            printf("-------[%d]--------\n",h);
            h++;
            APList apArray[3];
            map<int, APList >::iterator it = APMap.begin();
            int k = 0;
            for(int i=0;i < APMap.size() && k < M;i++,it++)
            {
                if(i == b[k])
                {
                    memcpy(&apArray[k], &it->second, sizeof(APList));
                    k++;
                }
            }
            mat X = CaculateByAPList(apArray);

            result += X;
        }
    }
}

int combineAPAndCacu()
{
    int sequ[3];
    combine(APMap.size(), 3, sequ,3);
    return 0;
}

int main(int argc, char ** argv)
{

    char apFile[256] = "./apinfo.116";
    char apFileENV[256] = "./ApEnv";

    if(argc > 1)
        strcpy(apFile, argv[1]);
    printf("file:%s\n",apFile);

    if(argc == 3)
        strcpy(apFileENV, argv[2]);
    printf("file:%s\n",apFileENV);

    if(CaculateAnByAPList(apFileENV))
    {
        printf("[%s-%d].ERROR: calu A and N fail!\n", __FILE__,__LINE__);
        return 0;
    }

    if(SetAPList(apFile))
    {
        printf("[%s-%d].ERROR: Set APlist faile!\n", __FILE__,__LINE__);
        return 0;
    }

    if(APMap.size() > 5)
    {
        map<int, APList >::iterator it = APMap.begin();
        for(int i =0 ;i <5; i++, it++);

        APMap.erase(it, APMap.end());
    }

    totalWeight = 0;
    combineAPAndCacu();
    //result.print("before deal result:");
    result = result / totalWeight;
    printf("\ncurr ENV A:%.2f N:%.2f\n",AP_A, AP_n);
    result.print("Last Result:");

    return 0;

}
