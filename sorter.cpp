//Code from https://mightynotes.wordpress.com/2019/03/21/when-python-is-not-enough-sorting-large-data-using-root-and-c/ 
void sorter(TString inpName, TString outName)
{
    //Open the previously created .root file
    TFile *f = new TFile(inpName);
    //Select the "Data" tree
    TTree *tree = (TTree*)f->Get("trimTree");
 
    //Create a list of sorting orders
    Int_t nentries = tree->GetEntries();
    Int_t *index = new Int_t[nentries];
 
    cout << "Done initalizing..." << endl;
 
    //Set how much data do we want to sort against (in this case, all data)
    tree->SetEstimate(nentries);
    //We are sorting against the "time" variable
    tree->Draw("energy_keV","","goff");
 
    cout << "Done loading list" << endl;
 
    //Perform the sorting
    TMath::Sort(nentries,tree->GetV1(),index);
    cout << "Done sorting" << endl;
     
 
    //Create a new file
    TFile f2(outName,"recreate");
    //Copy a empty tree so we have something to attach to
    TTree *tsorted = (TTree*)tree->CloneTree(0);
    //Setup caching and fetching parameters, making the process a lot faster
    tree->SetCacheSize(1000*1000*1000);//1G read cache
    tsorted->SetCacheSize(1000*1000*2000); //2G Write Cache
    tree->SetClusterPrefetch(true);
 
    auto start = chrono::high_resolution_clock::now();
    double last_time = 0;
    for (Int_t i=0;i<nentries;i++) {
        // Flush cached baskets to prevent OOM
        if((i+1)%0xfffff == 0)
            tree->DropBaskets();
        //Progress display
        if(i%8195 == 0) {
            auto now = chrono::high_resolution_clock::now();
            auto time_span = chrono::duration_cast<chrono::duration<double>>(now - start);
            double spend = time_span.count();
 
            double pct = 100.0*(double)i/nentries;
            printf("%c[2K\r", 27);
            cout << pct << "% "
                << i << "/" << nentries << ", "
                << "time: " << spend << "s. ETA: "
                << (size_t)((nentries-i)*(spend/i)) << "s"
                << ". delta = " << spend-last_time
                    << flush;
 
            last_time = spend;
        }
        //Fetch data from the source tree
        tree->GetEntry(index[i]);
        //And put it into the cloned tree
        tsorted->Fill();
    }
    cout << "\nDone sorting" << endl;
    tsorted->Write();
    delete [] index;
}
