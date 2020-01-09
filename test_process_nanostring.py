import unittest
from process_nanostring import parseRCC


class test_RCC(unittest.TestCase):
    """
    Tests for ouputs
    """
    def test_Sample_Level_Attributes(self):
        dfrcc, df_samp_attrib, df_lane_attrib = parseRCC(os.path.join(rcc_dir + "/", file))


        self.assertEqual()



if __name__ == '__main__':
    unittest.main()