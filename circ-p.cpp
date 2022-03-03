pids.push_back(fork());
srand48(time(0)+getpid());//seed random number generator


if(pids[pids.size()-1]){//if we're the child process

		{
			radii=list(2*sl*sl,.416666667,epsilonr);
	
		//makes one rectangular array, skewed by sqrtl(3)
			for(i=0;i<sl*sl;i++){
				if(i!=sl*sl/2+sl/2){
					//creates a normal particle
					x=((double)(i%sl));
					y=sqrtl(3.)*(i/sl);
					film.push_back(new particle(x,y,
						-h*(!dist(x,y,.4*sl,sl))
						,radii[radii.size()-1]
						,dist(x,y,.4*sl,sl),false));
				} else if(i==sl*sl/2+sl/2){
					//creates the particle depressed by the AFM
					film.push_back( new particle(
						(double)(i%sl)
						,sqrtl(3.)*(i/sl)
						,-h
						,radii[radii.size()-1]
						,true,true));
		  		} 
				radii.pop_back();
			}
			//makes another one, offset to complete the triangles.
			for(i=0;i<sl*sl;i++){
				//creates a normal particle
				x=((double)(i%sl))+.5;
				y=sqrtl(3.)*(i/sl)+sqrtl(3.)/2;
				film.push_back(new particle(x,y,
					-h*(!dist(x,y,.4*sl,sl),false)
					,radii[radii.size()-1]
					,dist(x,y,.4*sl,sl),false));
				radii.pop_back();
			}
	

			//add ligands
			for(i=0;i<sl;i++){
				for(j=0;j<sl;j++){
					//ligand in the j direction
					if(j!=sl-1){
						new ligand(film[sl*i+j],film[sl*i+j+1],prestr,1,2);
						new ligand(film[sl*sl+sl*i+j],film[sl*sl+sl*i+j+1],prestr,1,2);
					}
	
					//lower left to upper right ligand
					new ligand(film[sl*i+j],film[sl*sl+sl*i+j],prestr,1,2);
					if((j>0)){
						//lower right to upper left
						new ligand(film[sl*i+j],film[sl*sl-1+sl*i+j],prestr,1,2);
						//upper right to lower left
						if(i<sl-1){
							new ligand(film[sl+sl*i+j],film[sl*sl-1+sl*i+j],prestr,1,2);
						}
					}
					//upper left to lower right ligand
					if(i<sl-1){
						new ligand(film[sl+sl*i+j],film[sl*sl+sl*i+j],prestr,1,2);
					}
				}
			}
			todelete=deletelist(&film,parttoremove);//generate list
	if(fork()){
			deleter(&film,&todelete);//monte carlo element
			clean(&film);
			iterate(&film,h);//does the actual work.



			shapeu.open(("shape-circ-l-h"+Ttos(h)+"s"+Ttos(sl)+"p"+Ttos(prestr)+"#"+Ttos(k)+"a"+Ttos(afmr)+"q"+Ttos(getpid())+".txt").c_str(),ios::app);		
			forceu.open(("force-circ-l-h"+Ttos(h)+"s"+Ttos(sl)+"p"+Ttos(prestr)+"#"+Ttos(k)+"a"+Ttos(afmr)+"q"+Ttos(getpid())+".txt").c_str(),ios::app);		
			
			shapeu << film;
			forceu <<PRINTVAR <<", "<< multiforce(&film)<< endl;

			shapeu.close();
			forceu.close();
			for(i=0;i<film.size();i++){
				delete film[i];
			}
			film.clear();
		exit(0);
	}
		}
	//p=.1, from above
		{
			for(i=0;i<film.size();i++){
				delete film[i];
			}
			film.clear();
			//makes one rectangular array, skewed by sqrtl(3)
			for(i=0;i<sl*sl;i++){
				if(i!=sl*sl/2+sl/2){
					//creates a normal particle
					x=((double)(i%sl));
					y=sqrtl(3.)*(i/sl);
					filmb.push_back(new particle(x,y,0
						,radii[radii.size()-1]
						,dist(x,y,.4*sl,sl),false));
				} else if(i==sl*sl/2+sl/2){
					//creates the particle depressed by the AFM
					filmb.push_back( new particle(
						(double)(i%sl)
						,sqrtl(3.)*(i/sl)
						,-h
						,radii[radii.size()-1]
						,true,true));
			  	} 
				radii.pop_back();
			}
			//makes another one, offset to complete the triangles.
			for(i=0;i<sl*sl;i++){
				//creates a normal particle
				x=((double)(i%sl))+.5;
				y=sqrtl(3.)*(i/sl)+sqrtl(3.)/2;
				filmb.push_back(new particle(x,y,0
					,radii[radii.size()-1]
					,dist(x,y,.4*sl,sl),false));
				radii.pop_back();
			}

	
			//add ligands
			for(i=0;i<sl;i++){
				for(j=0;j<sl;j++){
					//ligand in the j direction
					if(j!=sl-1){
						new ligand(filmb[sl*i+j],filmb[sl*i+j+1],prestr,1,2);
						new ligand(filmb[sl*sl+sl*i+j],filmb[sl*sl+sl*i+j+1],prestr,1,2);
					}
	
					//lower left to upper right ligand
					new ligand(filmb[sl*i+j],filmb[sl*sl+sl*i+j],prestr,1,2);
					if((j>0)){
						//lower right to upper left
						new ligand(filmb[sl*i+j],filmb[sl*sl-1+sl*i+j],prestr,1,2);
						//upper right to lower left
						if(i<sl-1){
							new ligand(filmb[sl+sl*i+j],filmb[sl*sl-1+sl*i+j],prestr,1,2);
						}
					}
					//upper left to lower right ligand
					if(i<sl-1){
						new ligand(filmb[sl+sl*i+j],filmb[sl*sl+sl*i+j],prestr,1,2);
					}
				}
			}
			
			
			deleter(&filmb,&todelete);//monte carlo element
			clean(&filmb);
			iterate(&filmb,h);//does the actual work.
		
		
			shapeu.open(("shape-circ-u-h"+Ttos(h)+"s"+Ttos(sl)+"p"+Ttos(prestr)+"#"+Ttos(k)+"a"+Ttos(afmr)+"q"+Ttos(getpid())+".txt").c_str(),ios::app);		
			forceu.open(("force-circ-u-h"+Ttos(h)+"s"+Ttos(sl)+"p"+Ttos(prestr)+"#"+Ttos(k)+"a"+Ttos(afmr)+"q"+Ttos(getpid())+".txt").c_str(),ios::app);		
			
			shapeu << filmb;
			forceu <<PRINTVAR <<", "<<multiforce(&filmb)<< endl;

			shapeu.close();
			forceu.close();
			for(i=0;i<filmb.size();i++){
				delete filmb[i];
			}
			filmb.clear();
		}
		int * estatus=NULL;
		wait(estatus);
		exit(* estatus);
}

