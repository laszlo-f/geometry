	//p=0, from above
		{
			//makes one rectangular array, skewed by sqrt(3)
			for(i=0;i<sl*sl;i++){
				if(i!=sl*sl/2+sl/2){
					//creates a normal particle
					x=((long double)(i%sl));
					y=sqrtl(3.)*(i/sl);
					film.push_back(new particle(x,y,
						-h*(!dist(x,y,sl))
						,edge(i,sl)));
				} else if(i==sl*sl/2+sl/2){
					//creates the particle depressed by the AFM
					afm=new particle((long double)(i%sl),sqrtl(3.)*(i/sl),-h,true);
					film.push_back(afm);
			  	} 
			}
			//makes another one, offset to complete the triangles.
			for(i=0;i<sl*sl;i++){
				//creates a normal particle
				x=((long double)(i%sl))+.5;
				y=sqrtl(3.)*(i/sl)+sqrtl(3.)/2;
				film.push_back(new particle(x,y,
					-h*(!dist(x,y,sl))
					,edge(i,sl)));
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
			
			shapeu.open(("shape-s-l-p"+Ttos(prestr)+".txt").c_str(),ios::app);		
			forceu.open(("force-s-l-p"+Ttos(prestr)+".txt").c_str(),ios::app);		
			logforceu.open(("logforce-s-l-p"+Ttos(prestr)+".txt").c_str(),ios::app);		
			
			
			clean(&film);
			iterate(&film,h,afm);//does the actual work.
			shape << film;
			shapeu << film;
			force << h<< ", \t"<<afm->force()<< endl;
			forceu << h<< ", \t"<<afm->force()<< endl;
			logforce <<log(h)/log(10)<< ", \t"<<log(afm->force())/log(10)<< endl;
			logforceu <<log(h)/log(10)<< ", \t"<<log(afm->force())/log(10)<< endl;

			shapeu.close();
			forceu.close();
			logforceu.close();
			for(i=0;i<film.size();i++){
				delete film[i];
			}
			film.clear();
		}


