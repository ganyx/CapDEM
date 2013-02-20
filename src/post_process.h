#ifndef __POST_PROCESS_H__
#define __POST_PROCESS_H__

class Crun_post_process : public Crun //used to run parallel thread 
{
	public:
	std::vector <Cconfig> list_config; 		/**< List of configuration to be read, and stored*/
	int first;
	int last;
	string path_to_write;
	void init();
    void post_process();
  	virtual void run(){post_process(); };
};

#endif // __POST_PROCESS_H__
