//=========================================
// Function for reading in worldlines 
//=========================================
#include <builtin_types.h>
#include <stdio.h>
int getwl(float4 *worldlines, char *filename, int Nl, int Nppl, int Dim)
//Reads worldlines from filename and stores the results in worldlines.
//Worldline files may be prepared with wlGen.m in matlab
{
	int i,j,k;
	FILE *fp;
	printf("getting worldlines from %s\n",filename);
	if((fp = fopen(filename,"rb"))==NULL){
		printf("cannot open %s\n", *filename);
	}
	else{
		//loop over lines
		for(i=0;i<Nl;i++){
			//printf("getting worldline %d \n", i);
			//loop over points in line, x coordinate
			for(j=0;j<Nppl;j++){
				fscanf(fp,"%f ", &(worldlines[i*Nppl + j].x));
			}
			if(Dim>1){
				//loop over points in line, y coordinate
				for(j=0;j<Nppl;j++){
					fscanf(fp,"%f ",&(worldlines[i*Nppl + j].y));
				}
				if(Dim>2){
					//loop over points in line, z coordinate
					for(j=0;j<Nppl;j++){
						fscanf(fp,"%f ",&(worldlines[i*Nppl + j].z));
					}
				}
			}	
		}
	}
	if(fp) fclose(fp);  //Close file
	return 0;
}
