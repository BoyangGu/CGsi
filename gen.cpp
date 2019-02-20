//Created by: Adrian Diaz
//Last edited by: Adrian Diaz
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#define PI 3.14159265
using namespace std;

int main()
{
	string line;
	int aiuc; //number of atoms in the unit cell
	int aiuc2;
	double cell[3][3]; //unit cell basis vectors and nanomesh basis vectors
	double cell2[3][3];
	double px, py, pz, pxn, pyn, pzn;
	double** atom_positions; //Will store the positions of the atoms in the unit cell
	int repeat[3]; //repeats of the unit cell in each dimension,
	double scale = 10;           //number of unit cells in one element
    int elements{0};//number of elements in the specific region
    int numele{0};//restore the node position of element 
    double origin[3];
    double originele[3];
    double origincell[3];
	//compared to primitive unit cell, assumes orthogonal primitive unit cell
	int i, j, k, l, m, mm, mmm, atomn; //counter variables
	double n1x, n1y, n1z, n2x, n2y, n2z, n3x, n3y, n3z, n4x, n4y, n4z
		, n5x, n5y, n5z, n6x, n6y, n6z, n7x, n7y, n7z, n8x, n8y, n8z;
	int *atom_types;
	double *atom_charges;
//    double xlo{0},xhi{0},ylo{0},yhi{0},zlo{0},zhi{0};

	ofstream myfile("Si_Model.txt"); //output filestream object for file output
	ifstream infile("Basis.si"); //input filestream object for file reading
	if (infile.is_open()) {
		infile >> line;
		cout << line << " ";
		infile >> line;
		cout << line << " ";
		infile >> line;
		cout << line << "\n";
		infile >> aiuc;
		aiuc2 = aiuc;
		cout << aiuc;
		
		for (i = 0; i<3; i++) {
			cout << "\n";
			for (j = 0; j<3; j++) {
				infile >> cell[i][j];
				cout << cell[i][j] << " ";
				cell2[i][j] = cell[i][j];
				cell[i][j] = scale * cell[i][j];
			}
		}
		cout << "\n";

		infile >> repeat[0];  

		infile >> repeat[1];

		infile >> repeat[2];


		cout << repeat[0] << " " << repeat[1] << " " << repeat[2] << "\n";
		atom_positions = new double*[aiuc2];
		atom_types = new int[aiuc2];
		for (i = 0; i<aiuc2; i++) {
			atom_positions[i] = new double[3];
		}
		for (i = 0; i<aiuc2; i++) {
			infile >> atom_types[i];
			cout << line << " ";
			infile >> atom_positions[i][0];
			infile >> atom_positions[i][1];
			infile >> atom_positions[i][2];
			cout << atom_positions[i][0] << " " << atom_positions[i][1] << " " << atom_positions[i][2] << "\n";
		}

	}


	else {
		cout << "unable to open file";
	}
//    origincell[0]=(cell2[0][0]+cell2[1][0]+cell2[2][0])/2;
//    origincell[1]=(cell2[0][1]+cell2[1][1]+cell2[2][2])/2;
//    origincell[2]=(cell2[0][2]+cell2[1][2]+cell2[2][2])/2;
    originele[0]=(cell[0][0]+cell[1][0]+cell[2][0])/2;
    originele[1]=(cell[0][1]+cell[1][1]+cell[2][1])/2;
    originele[2]=(cell[0][2]+cell[1][2]+cell[2][2])/2;
    origin[0]=(cell[0][0]*repeat[0]+cell[1][0]*repeat[1]+cell[2][0]*repeat[2])/2;   //determine the center of the model
    origin[1]=(cell[0][1]*repeat[0]+cell[1][1]*repeat[1]+cell[2][1]*repeat[2])/2;
    origin[2]=(cell[0][2]*repeat[0]+cell[1][2]*repeat[1]+cell[2][2]*repeat[2])/2;
//    origin[0]=0;
//    origin[1]=0;
//    origin[2]=0;
    cout<<"repeat in x: "<<repeat[0]<<endl;
    cout<<"repeat in y: "<<repeat[1]<<endl;
    cout<<"repeat in z: "<<repeat[2]<<endl;

    cout<<"origin[0]: "<<origin[0]<<endl;
    cout<<"origin[1]: "<<origin[1]<<endl;
    cout<<"origin[2]: "<<origin[2]<<endl;
    
        		for (m = 0; m<repeat[0] ; m++) {                    //in the ith element
                 for (mm = 0; mm<repeat[1] ; mm++) {
                  for (mmm = 0; mmm<repeat[2] ; mmm++){
										pxn = cell[0][0] * (m-1) + cell[1][0] * (mm-1) + cell[2][0] * (mmm-1)-origin[0]+ originele[0];             //center of element
										pyn = cell[0][1] * (m-1) + cell[1][1] * (mm-1) + cell[2][1] * (mmm-1)-origin[1]+originele[1] ;
										pzn = cell[0][2] * (m-1) + cell[1][2] * (mm-1) + cell[2][2] * (mmm-1)-origin[2]+originele[2] ;
                                            if(pxn*pxn>=400||pzn<=500){     
//                                              if((pxn<=-40||pxn>=-50)&&pzn*pzn<=52900){											
                                                elements+=1;
											  }
//                                            }
                 }
                }
               }
         cout<<"number of elements: "<<elements;      
	if (myfile.is_open()) {
		myfile << fixed << setprecision(8);
		myfile << "Si \n \n";
		myfile << elements << " " << "CAC_elements \n";
		myfile << "1 atom types \n \n";
//		myfile << (cell[0][0]*0   -0.5) <<"  " << cell[0][0] * repeat[0] + 0.5 << "  xlo xhi \n";
//		myfile << (cell[1][1]*0   -0.5) <<"  " << cell[1][1] * repeat[1] + 0.5 << "  ylo yhi \n";
//		myfile << (cell[2][2]*0  -0.5) <<"  " << cell[2][2] * repeat[2] + 0.5 << "  zlo zhi \n";
		myfile << -1550 <<"  " << 1550 << "  xlo xhi \n";
		myfile << -1550 <<"  " << 1550 << "  ylo yhi \n";
		myfile << -1550 <<"  " << 1550 << "  zlo zhi \n";
		myfile << "Masses \n \n";
		myfile << "1          28.0855\n";
		myfile << "CAC_Elements \n \n";
        atomn=0;
		for (m = 0; m<repeat[0] ; m++) {                    //in the ith element
			for (mm = 0; mm<repeat[1] ; mm++) {
				for (mmm = 0; mmm<repeat[2] ; mmm++) {
										pxn = cell[0][0] * (m-1) + cell[1][0] * (mm-1) + cell[2][0] * (mmm-1)-origin[0]+originele[0] ;
										pyn = cell[0][1] * (m-1) + cell[1][1] * (mm-1) + cell[2][1] * (mmm-1)-origin[1]+originele[1] ;
										pzn = cell[0][2] * (m-1) + cell[1][2] * (mm-1) + cell[2][2] * (mmm-1)-origin[2] +originele[2];
                                            if(pxn*pxn>=400||pzn<=500){     
//                                               if((pxn<=-40||pxn>=-50)&&pzn*pzn<=52900){											
											
                                        myfile << std::left << 
                                            setw(15) << atomn+1 << "Eight_Node  " << setw(15) << aiuc2 << setw(15) << int(scale) << setw(15) << int(scale) <<setw(15) << int(scale) << "\n" ;                                                                          
                                        numele+=1;       									
                                    for (l = 0; l < aiuc2; l++) {
										px = cell[0][0] * m + cell[1][0] * mm + cell[2][0] * mmm + atom_positions[l][0]-origin[0];
										py = cell[0][1] * m + cell[1][1] * mm + cell[2][1] * mmm + atom_positions[l][1]-origin[1];
										pz = cell[0][2] * m + cell[1][2] * mm + cell[2][2] * mmm + atom_positions[l][2]-origin[2];


										n1x = px - (scale) / 2 * cell2[0][0] - (scale) / 2 * cell2[1][0] - (scale) / 2 * cell2[2][0];
										n1y = py - (scale) / 2 * cell2[0][1] - (scale) / 2 * cell2[1][1] - (scale) / 2 * cell2[2][1];
										n1z = pz - (scale) / 2 * cell2[0][2] - (scale) / 2 * cell2[1][2] - (scale) / 2 * cell2[2][2];

										n2x = px + (scale) / 2 * cell2[0][0] - (scale) / 2 * cell2[1][0] - (scale) / 2 * cell2[2][0];
										n2y = py + (scale) / 2 * cell2[0][1] - (scale) / 2 * cell2[1][1] - (scale) / 2 * cell2[2][1];
										n2z = pz + (scale) / 2 * cell2[0][2] - (scale) / 2 * cell2[1][2] - (scale) / 2 * cell2[2][2];
                                        
										n3x = px + (scale) / 2 * cell2[0][0] + (scale) / 2 * cell2[1][0] - (scale) / 2 * cell2[2][0];
										n3y = py + (scale) / 2 * cell2[0][1] + (scale) / 2 * cell2[1][1] - (scale) / 2 * cell2[2][1];
										n3z = pz + (scale) / 2 * cell2[0][2] + (scale) / 2 * cell2[1][2] - (scale) / 2 * cell2[2][2];

										n4x = px - (scale) / 2 * cell2[0][0] + (scale) / 2 * cell2[1][0] - (scale) / 2 * cell2[2][0];
										n4y = py - (scale) / 2 * cell2[0][1] + (scale) / 2 * cell2[1][1] - (scale) / 2 * cell2[2][1];
										n4z = pz - (scale) / 2 * cell2[0][2] + (scale) / 2 * cell2[1][2] - (scale) / 2 * cell2[2][2];

										n5x = px - (scale) / 2 * cell2[0][0] - (scale) / 2 * cell2[1][0] + (scale) / 2 * cell2[2][0];
										n5y = py - (scale) / 2 * cell2[0][1] - (scale) / 2 * cell2[1][1] + (scale) / 2 * cell2[2][1];
										n5z = pz - (scale) / 2 * cell2[0][2] - (scale) / 2 * cell2[1][2] + (scale) / 2 * cell2[2][2];
                                        
										n6x = px + (scale) / 2 * cell2[0][0] - (scale) / 2 * cell2[1][0] + (scale) / 2 * cell2[2][0];
										n6y = py + (scale) / 2 * cell2[0][1] - (scale) / 2 * cell2[1][1] + (scale) / 2 * cell2[2][1];
										n6z = pz + (scale) / 2 * cell2[0][2] - (scale) / 2 * cell2[1][2] + (scale) / 2 * cell2[2][2];

										n7x = px + (scale) / 2 * cell2[0][0] + (scale) / 2 * cell2[1][0] + (scale) / 2 * cell2[2][0];
										n7y = py + (scale) / 2 * cell2[0][1] + (scale) / 2 * cell2[1][1] + (scale) / 2 * cell2[2][1];
										n7z = pz + (scale) / 2 * cell2[0][2] + (scale) / 2 * cell2[1][2] + (scale) / 2 * cell2[2][2];

										n8x = px - (scale) / 2 * cell2[0][0] + (scale) / 2 * cell2[1][0] + (scale) / 2 * cell2[2][0];
										n8y = py - (scale) / 2 * cell2[0][1] + (scale) / 2 * cell2[1][1] + (scale) / 2 * cell2[2][1];
										n8z = pz - (scale) / 2 * cell2[0][2] + (scale) / 2 * cell2[1][2] + (scale) / 2 * cell2[2][2];
										myfile << std::left <<
											setw(15) << 1 << setw(15) << l + 1 << setw(15) << atom_types[l] << setw(25) << n1x << setw(25) << n1y << setw(25) << n1z << "\n" <<
											setw(15) << 2 << setw(15) << l + 1 << setw(15) << atom_types[l] << setw(25) << n2x << setw(25) << n2y << setw(25) << n2z << "\n" <<
											setw(15) << 3 << setw(15) << l + 1 << setw(15) << atom_types[l] << setw(25) << n3x << setw(25) << n3y << setw(25) << n3z << "\n" <<
											setw(15) << 4 << setw(15) << l + 1 << setw(15) << atom_types[l] << setw(25) << n4x << setw(25) << n4y << setw(25) << n4z << "\n" <<
											setw(15) << 5 << setw(15) << l + 1 << setw(15) << atom_types[l] << setw(25) << n5x << setw(25) << n5y << setw(25) << n5z << "\n" <<
											setw(15) << 6 << setw(15) << l + 1 << setw(15) << atom_types[l] << setw(25) << n6x << setw(25) << n6y << setw(25) << n6z << "\n" <<
											setw(15) << 7 << setw(15) << l + 1 << setw(15) << atom_types[l] << setw(25) << n7x << setw(25) << n7y << setw(25) << n7z << "\n" <<
											setw(15) << 8 << setw(15) << l + 1 << setw(15) << atom_types[l] << setw(25) << n8x << setw(25) << n8y << setw(25) << n8z << "\n";
                                    }
                                            atomn+=1;                                    
									}
//									}
								}
							}
						}
		cout << "done\n";					
}
	else {
		cout << "Unable to open file";
	}
	/*
	if (myfile2.is_open()) {
		myfile2 << "Hello " << "        \n";//arbitrary lineholder for the rewrite of the number of atoms later on
		myfile2 << "Atoms. Timestep: 0 \n";

		atomn = 1;
		for (m = 0; m<repeat[0] / Nx; m++) {
			cout << m << "\n";
			for (mm = 0; mm<repeat[1] / Ny; mm++) {
				for (mmm = 0; mmm<repeat[2] / Nz; mmm++) {
					for (i = 0; i<Nx; i++) {
						for (j = 0; j<Ny; j++) {
							for (k = 0; k<Nz; k++) {
								for (l = 0; l<aiuc2; l++) {
									px = meshcell[0][0] * m + meshcell[1][0] * mm + meshcell[2][0] * mmm +
										cell[0][0] * i + cell[1][0] * j + cell[2][0] * k + atom_positions[l][0];
									py = meshcell[0][1] * m + meshcell[1][1] * mm + meshcell[2][1] * mmm +
										cell[0][1] * i + cell[1][1] * j + cell[2][1] * k + atom_positions[l][1];
									pz = meshcell[0][2] * m + meshcell[1][2] * mm + meshcell[2][2] * mmm +
										cell[0][2] * i + cell[1][2] * j + cell[2][2] * k + atom_positions[l][2];





									
										myfile2 << std::left << "1  " << setw(10) << px << setw(10) << py << setw(10) << pz << "\n";
										myfile3 << std::left << setw(10) << atomn << setw(10) << l + 1 << setw(10) << px << setw(10) << py << setw(10) << pz << "\n";
										atomn = atomn + 1;
									





								}
							}
						}
					}
				}
			}
		}
		myfile2.seekp(0, std::ios::beg);
		myfile2 << atomn - 1;
		cout << atomn - 1 << "\n";
	}
	else {
		cout << "Unable to open file";
	}
	*/
	for (i = 0; i<3; i++) {

		for (j = 0; j<3; j++) {
			cout << cell[i][j] << " ";
		}
		cout << "\n";
	}
	return 0;
}
