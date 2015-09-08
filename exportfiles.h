


//a list which txt-files to produce
    ofstream mikadofile("mikado.txt"); 
    ofstream nodefile("nodes.txt");
    ofstream springfile("springs.txt");
    ofstream XYfile("conjpoints.txt");
    ofstream shearcoordinates("shearcoordinates.txt");
    ofstream shearenergy("shearenergy.txt");
    ofstream cluster("clusters.txt");
    ofstream clusterdata("clusterdata.txt", ios_base::app | ios_base::out);
    ofstream anglefile("angles.txt");
    char s[80];
    ofstream clusterdistribution(s);