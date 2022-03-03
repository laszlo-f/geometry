pids.push_back(fork());
srand48(time(0)+getpid());//seed random number generator
srand(time(0)+getpid());//seed random number generator

if(pids[pids.size()-1]){//if we're the child process
	radii=list(2*sl*sl,rscale(naturalradius,naturalradius,naturalgap),epsilonr);//make a list of radii
	radiib=radii;//duplicate of list must be created because first list is destroyed before fork()
	//I should probably be using iterators with these lists instead of pop_back() to avoid this issue
	
	//makes one rectangular array, skewed by sqrtl(3)
	for(i=0;i<sl*sl;i++){
		//creates a normal particle
		x=((double)(i%sl));
		y=sqrtl(3.)*(i/sl);
		film.push_back(new particle( x ,y
			,-h*(!dist(x,y,.4*sl,sl))
			,radii[radii.size()-1]
			, dist(x,y,.4*sl,sl)||!dist(x-offset,y,afmr*sl,sl)
			,!dist(x-offset,y,afmr*sl,sl)));
		radii.pop_back();
	}
	//makes another one, offset to complete the triangles.
	for(i=0;i<sl*sl;i++){
		//creates a normal particle
		x=((double)(i%sl))+.5;
		y=sqrtl(3.)*(i/sl)+sqrtl(3.)/2;
		film.push_back(new particle( x ,y
			,-h*(!dist(x,y,.4*sl,sl))
			,radii[radii.size()-1]
			, dist(x,y,.4*sl,sl)||!dist(x-offset,y,afmr*sl,sl)
			,!dist(x-offset,y,afmr*sl,sl)));
		radii.pop_back();
	}
	spconsts=list(3*film.size(),1,epsilon);//make a list of spring constants
	unstrech=list(3*film.size(),unstrecher(prestr,naturalradius,naturalradius,naturalgap),epsilonu);//make a list of unstretched lengths
	todelete=deletelist(&film,DELETEARGS);//generate list of particles to delete
	//fork() must come after deletelist because its output is used in both branches
	//deletelist needs the constructed film because it only wants to delete unifixed particles.



///////////////////////////////////////////////////////


	//in this branch, we move particles up from below using the already created vector
	if(fork()){
		//add ligands
		for(i=0;i<sl;i++){
			for(j=0;j<sl;j++){
				//ligand in the j direction
				if(j!=sl-1){
					new ligand(film[sl*i+j],film[sl*i+j+1],unstrech[unstrech.size()-1],spconsts[spconsts.size()-1],2);
					spconsts.pop_back();
					unstrech.pop_back();
					new ligand(film[sl*sl+sl*i+j],film[sl*sl+sl*i+j+1],unstrech[unstrech.size()-1],spconsts[spconsts.size()-1],2);
					spconsts.pop_back();
					unstrech.pop_back();
				}
	
				//lower left to upper right ligand
				new ligand(film[sl*i+j],film[sl*sl+sl*i+j],unstrech[unstrech.size()-1],spconsts[spconsts.size()-1],2);
					spconsts.pop_back();
					unstrech.pop_back();
				if((j>0)){
					//lower right to upper left
					new ligand(film[sl*i+j],film[sl*sl-1+sl*i+j],unstrech[unstrech.size()-1],spconsts[spconsts.size()-1],2);
					spconsts.pop_back();
					unstrech.pop_back();
					//upper right to lower left
					if(i<sl-1){
						new ligand(film[sl+sl*i+j],film[sl*sl-1+sl*i+j],unstrech[unstrech.size()-1],spconsts[spconsts.size()-1],2);
						spconsts.pop_back();
						unstrech.pop_back();
					}
				}
				//upper left to lower right ligand
				if(i<sl-1){
					new ligand(film[sl+sl*i+j],film[sl*sl+sl*i+j],unstrech[unstrech.size()-1],spconsts[spconsts.size()-1],2);
					spconsts.pop_back();
					unstrech.pop_back();
				}
			}
		}
		
		
		deleter(&film,&todelete);//monte carlo element
					//removes random particles
					//listed in todelete
		clean(&film);//remove unneeded data that will slow down iterate()
		iterate(&film,h);//does the actual work.

		//open files for output
		shapeu.open((
			"shape-circ-l-h"+Ttos(h)
			+"s"+Ttos(sl)
			+"p"+Ttos(prestr)
			+"ek"+Ttos(epsilon)
			+"eu"+Ttos(epsilonu)
			+"er"+Ttos(epsilonr)
			+"ed"+Ttos(epsilond)
			+"#"+Ttos(k)
			+"a"+Ttos(afmr)
			+"o"+Ttos(offset)
			+"n"+Ttos(naturalradius)
			+"c"+Ttos(curve(&film))
			+"f"+Ttos(multiforce(&film))
			+"q"+Ttos(getpid())
			+".txt"
		).c_str(),ios::app);		
		forceu.open((
			"force-circ-l-h"+Ttos(h)
			+"s"+Ttos(sl)
			+"p"+Ttos(prestr)
			+"ek"+Ttos(epsilon)
			+"eu"+Ttos(epsilonu)
			+"er"+Ttos(epsilonr)
			+"ed"+Ttos(epsilond)
			+"#"+Ttos(k)
			+"a"+Ttos(afmr)
			+"o"+Ttos(offset)
			+"n"+Ttos(naturalradius)
			+"c"+Ttos(curve(&film))
			+"f"+Ttos(multiforce(&film))
			+"q"+Ttos(getpid())
			+".txt"
		).c_str(),ios::app);		
		
		//do outputing
		shapeu << film;
		lowerupper="lower";
		forceu <<h<<", "<<sl<<", "<<prestr<<", "<<epsilon<<", "<<epsilonu<<", "<<epsilonr<<", "<<epsilond<<", "<<k<<", "<<afmr<<", "<<offset<<", "<<naturalradius<<", "<<lowerupper<<", "<<getpid()<<", "<<curve(&film)<<", " << multiforce(&film)<< endl;

		//close files
		shapeu.close();
		forceu.close();

		//free memory
		for(i=0;i<film.size();i++){
			delete film[i];
		}
		film.clear();
		exit(0);//quit with success!
	}

//////////////////////////////////////////////////////////////
	//In this branch, we move particles down from above
	//
	for(i=0;i<film.size();i++){//free memory
		delete film[i];
	}
	film.clear();

	//makes one rectangular array, skewed by sqrtl(3)
	for(i=0;i<sl*sl;i++){
		//creates a normal particle
		x=((double)(i%sl));
		y=sqrtl(3.)*(i/sl);
		filmb.push_back(new particle( x ,y
			,-h*(!dist(x-offset,y,afmr*sl,sl))
			,radiib[radiib.size()-1]
			, dist(x,y,.4*sl,sl)||!dist(x-offset,y,afmr*sl,sl)
			,!dist(x-offset,y,afmr*sl,sl)));
		radiib.pop_back();
	}
	//makes another one, offset to complete the triangles.
	for(i=0;i<sl*sl;i++){
		//creates a normal particle
		x=((double)(i%sl))+.5;
		y=sqrtl(3.)*(i/sl)+sqrtl(3.)/2;
		filmb.push_back(new particle( x ,y
			,-h*(!dist(x-offset,y,afmr*sl,sl))
			,radiib[radiib.size()-1]
			, dist(x,y,.4*sl,sl)||!dist(x-offset,y,afmr*sl,sl)
			,!dist(x-offset,y,afmr*sl,sl)));
		radiib.pop_back();
	}
	
	//add ligands
	for(i=0;i<sl;i++){
		for(j=0;j<sl;j++){
			//ligand in the j direction
			if(j!=sl-1){
				new ligand(filmb[sl*i+j],filmb[sl*i+j+1],unstrech[unstrech.size()-1],spconsts[spconsts.size()-1],2);
				spconsts.pop_back();
				unstrech.pop_back();
				new ligand(filmb[sl*sl+sl*i+j],filmb[sl*sl+sl*i+j+1],unstrech[unstrech.size()-1],spconsts[spconsts.size()-1],2);
				spconsts.pop_back();
				unstrech.pop_back();
			}
	
			//lower left to upper right ligand
			new ligand(filmb[sl*i+j],filmb[sl*sl+sl*i+j],unstrech[unstrech.size()-1],spconsts[spconsts.size()-1],2);
			spconsts.pop_back();
						unstrech.pop_back();
			if((j>0)){
				//lower right to upper left
				new ligand(filmb[sl*i+j],filmb[sl*sl-1+sl*i+j],unstrech[unstrech.size()-1],spconsts[spconsts.size()-1],2);
				spconsts.pop_back();
				unstrech.pop_back();
				//upper right to lower left
				if(i<sl-1){
					new ligand(filmb[sl+sl*i+j],filmb[sl*sl-1+sl*i+j],unstrech[unstrech.size()-1],spconsts[spconsts.size()-1],2);
					spconsts.pop_back();
					unstrech.pop_back();
				}
			}
			//upper left to lower right ligand
			if(i<sl-1){
				new ligand(filmb[sl+sl*i+j],filmb[sl*sl+sl*i+j],unstrech[unstrech.size()-1],spconsts[spconsts.size()-1],2);
				spconsts.pop_back();
				unstrech.pop_back();
			}
		}
	}
			
	deleter(&filmb,&todelete);//monte carlo element
	clean(&filmb);
	iterate(&filmb,h);//does the actual work.
		
	shapeu.open((
		"shape-circ-u-h"+Ttos(h)
		+"s"+Ttos(sl)
		+"p"+Ttos(prestr)
		+"ek"+Ttos(epsilon)
		+"eu"+Ttos(epsilonu)
		+"er"+Ttos(epsilonr)
		+"ed"+Ttos(epsilond)
		+"#"+Ttos(k)
		+"a"+Ttos(afmr)
		+"o"+Ttos(offset)
		+"n"+Ttos(naturalradius)
		+"c"+Ttos(curve(&filmb))
		+"f"+Ttos(multiforce(&filmb))
		+"q"+Ttos(getpid())
		+".txt"
	).c_str(),ios::app);		
	forceu.open((
		"force-circ-u-h"+Ttos(h)
		+"s"+Ttos(sl)
		+"p"+Ttos(prestr)
		+"ek"+Ttos(epsilon)
		+"eu"+Ttos(epsilonu)
		+"er"+Ttos(epsilonr)
		+"ed"+Ttos(epsilond)
		+"#"+Ttos(k)
		+"a"+Ttos(afmr)
		+"o"+Ttos(offset)
		+"n"+Ttos(naturalradius)
		+"c"+Ttos(curve(&filmb))
		+"f"+Ttos(multiforce(&filmb))
		+"q"+Ttos(getpid())
		+".txt"
	).c_str(),ios::app);		
			
	shapeu << filmb;
	lowerupper="upper";
	forceu <<h<<", "<<sl<<", "<<prestr<<", "<<epsilon<<", "<<epsilonu<<", "<<epsilonr<<", "<<epsilond<<", "<<k<<", "<<afmr<<", "<<offset<<", "<<naturalradius<<", "<<lowerupper<<", "<<getpid()<<", "<<curve(&filmb)<<", " << multiforce(&filmb)<< endl;

	shapeu.close();
	forceu.close();
	for(i=0;i<filmb.size();i++){
		delete filmb[i];
	}
	filmb.clear();
	int * estatus=NULL;
	wait(estatus);//wait for child to exit
	exit(* estatus);//tell parent we're done
}

