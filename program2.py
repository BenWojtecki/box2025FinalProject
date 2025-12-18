from collections import Counter
import pandas as pan

def naive(maws, threshold):
    """
    maws given as a list
    """
    mawsDF = pan.DataFrame(Counter(maws).items(),columns=("maw","occ"))

    return mawsDF[mawsDF["occ"]>threshold].sort_values(by=["occ","maw"],ascending=[False,True])


print(naive(["a","a","a","a","b","b","b","c","z","z","z","z"],2))