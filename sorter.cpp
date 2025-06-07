void sorter(TString inpName, TString outName)
{
    // Open the input file and get the tree
    TFile *fin = new TFile(inpName, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Error opening input file: " << inpName << std::endl;
        return;
    }

    TTree *tree = (TTree*)fin->Get("trimTree");
    if (!tree) {
        std::cerr << "Error: tree 'trimTree' not found in file." << std::endl;
        return;
    }

    std::cout << "Building index on energy_eV..." << std::endl;
    tree->BuildIndex("energy_eV");  // use index for sorting

    // Open output file
    TFile fout(outName, "RECREATE");
    if (fout.IsZombie()) {
        std::cerr << "Error creating output file: " << outName << std::endl;
        return;
    }

    std::cout << "Copying sorted tree..." << std::endl;
    TTree *sorted = tree->CopyTree("");  // uses the index
    sorted->Write("trimTree", TObject::kOverwrite);

    std::cout << "Done." << std::endl;

    delete fin;
}
