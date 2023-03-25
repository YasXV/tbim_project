#include<iostream>
#include<cmath>
#include<fstream>
#include<string>
#include<vector>
#include<sstream>
#include<typeinfo>
#include<regex>
#include<algorithm>
#include<iterator>
#include <sys/stat.h>
#include <sys/types.h>
#include <iomanip>
#include <ctime>
using namespace std;

#ifndef __nano_h__
#define __nano_h__

namespace nano{
  //test of opening the file-----------------
  bool open_file(string filename);

  //split a string by spaces-----------------
  vector<string> split(const string com);

  //construction of all parameters (Needed for starting calcul)-----------
  void initialize(string filename);

  //randomizer-------------------------------
  int randomf(int range);

  //create the name of the file of the output:
  string name_file(int mnomi);

  //Function to calculate distance-----------
  double distance(float xij,float yij, float zij);

  //Function to write some parameters to plot courb in the future :
  void write_courbe(string filename);

  //Classe of nanoparticle-------------------
  class Nanoparticle{
  private:
    float apar; //=-0.0135;
    float ecoh; //=0.55/12.;
    float x,y,z;
    string kind;
    string type_vois;
  public:
    //creator---------------------------------
    Nanoparticle();

    //setting the parameters in the maille (kind,x,y,z) :
    void set_param(string line);

    //return position in maille :--------------
    vector<float> position();

    //change the kind of the nanoparticule-----
    void change_type(string new_type);

    //return the kind of the particle:---------
    const string getKind();

    //
    const string getType_vois();

    //Destructor------------------------------
    ~Nanoparticle();
  };

  //Classe Maille-----------------------------
  class Maille{
  private:
    vector<Nanoparticle> maille;
    vector<int> nvois;
    vector<vector<int>> ivois;
    float N_atoms;
    vector<float> energy_atoms;
  public:
    //creation de la structure ----------------
    Maille(string filename);

    //create vectors ivois and nvois :---------
    void voisin();

    //return of ivois
    const vector<vector<int>> getIVois();

    //change the kind of particle i :----------
    void changeParticle(int site,string new_kind);

    //Function to write all parameters in a file :
    void write_parameters(string filename,int pas);

    //return the entire maille-----------------
    const vector<Nanoparticle> getMaille();

    //return a nanoparticle of the maille------
    const Nanoparticle getParticle(int position);

    //return of nvois
    const vector<int> getNVois();

    //return the number of atom in the maille--
    const float getNumberOfAtoms();

    //return the kind of the particle i :------
    const string getParticleKind(int sitePosition);

    //return the type of particuke(its placement in the maille) depending on the number of neighbours it has
    const string getParticleType(int sitePosition);

    //surcharge de l'operateur = :-------------
    Maille operator=(const Maille& source);

    //Destructor-------------------------------
    ~Maille();
  };

//calulate the energy difference
  float diffenergy(Maille& maille,int site,vector<float>& ener_0) ;

  //change an atom of the maille in the impurity :
  int mc_exchange(Maille& maille,int ipas);

  //Monte Carlo of a maille (maybe to put in class maille)---------
  vector<float> Monte_Carlo(Maille& maille);

  void DoMonteCarlo(Maille& maille);

  //count the number of atoms depending on their types and stock it into a vector of size 5(5 different types)
  void count_type_atoms(vector<int> &numb_type, Maille maille);
  //all elements of a vector becomes 0
  void clear_vect(vector<float> vect);

  //write in a file one concentration
  void writeConcen(string path,float concen);

  //creates all the concentrations files where we write the concentration of all types of particule
  void create_concen_files(string directory);

  //write the concentration for all type of particules
  void write_all(string directory,vector<float> c2meanvect);

  //determines the concentration for each type of particule(surface, heart..)
  void concen_neighbor(Maille maille,vector<float> &c2meanvect);

}
#endif
