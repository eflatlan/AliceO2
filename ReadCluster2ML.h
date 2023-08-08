class ReadCluster2ML
{
 public:
  ReadCluster2ML() = default;
  //  : mReadFile(readFile) {}
  ~ReadCluster2ML() override = default;

  void init(framework::InitContext& ic) final;

  void run(framework::ProcessingContext& pc) final;
  // void endOfStream(framework::EndOfStreamContext& ec) override;

 private:
  void initFileIn(const std::string& filename);

  long mClustersReceived;
  ExecutionTimer mExTimer;

  std::unique_ptr<TFile> mFile; // root file with Clusters
  std::unique_ptr<TTree> mTree; // tree inside the file
  std::vector<o2::hmpid::Trigger> mClusterTriggersFromFile, *mClusterTriggersFromFilePtr = &mClusterTriggersFromFile;
  std::vector<o2::hmpid::Cluster> mClustersFromFile, *mClustersFromFilePtr = &mClustersFromFile;

  unsigned long mNumberOfEntries = 0; // number of entries from TTree
  unsigned long mCurrentEntry = 0;    // index of current entry

  // void strToFloatsSplit(std::string s, std::string delimiter, float* res, int maxElem = 7);
};