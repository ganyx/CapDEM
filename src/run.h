class Cevent{
  public:
//    double period, next;
    int period, next;
    bool should_do(double t);
    
};

/*! \class Crun 
*  \brief This class is the master class in the Soft-Dynamics program. It is used to know where to read an initial configuration, where and when to save it, what to do (simulate the grains dynamics or post-process) and for how long.
*/

class Crun //: public QThread //used to run parallel thread 
{
  public:
    class Cin_out   where_read;      /**<Set of files where to read the configuration. */
    class Cin_out   where_save;     /**<Set of files where to save the configuration.*/

    class Cconfig config;          /**<Configuration (particles, contacts, cell, physical parameters).*/

    class Cevent screen;           /**<Define when should print of screen some simulation details.*/
    class Cevent save;             /**<Define when should save on file the configuration.*/
  //  class Cpost_process post_process;
  
    double t;                     /**< Current time.*/
    double tstart;                /**< Initial time at the begening of the simulation.*/
    double tend;                  /**< Time when the sinmulation will end.*/
    double dt;                    /**< Time step.*/
   	 
   	 
   // void init_create();
    void init_evolve();
    void evolve();                /**< Iterate the config between tstart and tend with dt step.*/ 
   	
};



class Crun_create : public Crun //used to run parallel thread 
{
  virtual void run(){ };
};

class Crun_evolve : public Crun //used to run parallel thread 
{
  virtual void run(){evolve();};
};













