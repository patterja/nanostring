import unittest
from process_nanostring import parseRCC
class test_RCC(unittest.TestCase):
    """
    Tests for ouputs
    """
    def test_Sample_Level_Attributes(self):
        dfrcc, df_samp_attrib, df_lane_attrib = parseRCC(os.path.join(rcc_dir + "/", file))


        self.assertEqual()
    for v in range(len(samp_attrib_dict.values())-1):
        if samp_attrib_dict.values()[v].equals(samp_attrib_dict.values()[v+1]) != True:
            print("Samples are not from one batch. Sample Attributes differ")
            print(samp_attrib_dict.keys()[v])
        else:
            print("Samples are from one batch. OK")


    def (self):
        self.assertEqual(sum([1, 2, 3]), 6, "Should be 6")

    def test_sum_tuple(self):
        self.assertEqual(sum((1, 2, 2)), 6, "Should be 6")

