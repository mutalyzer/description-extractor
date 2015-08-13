"""
Tests for the mutalyzer.describe module.
"""


from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from extractor import describe_protein
from extractor.util import str


class TestDescribe:
    """
    Test the mutalyzer.describe module.
    """
    def setup(self):
        self.reference = ('MAVLWRLSAVCGALGGRALLLRTPVVRPAHISAFLQDRPIPEWCGVQHI' +
            'HLSPSHHSGSKAASLHWTSERVVSVLLLGLLPAAYLNPCSAMDYSLAAALTLHGHWGLGQVVT' +
            'DYVHGDALQKAAKAGLLALSALTFAGLCYFNYHDVGICKAVAMLWKL*')
        self.sample = [
            ('MAVLWRLVCGALGGRALLLRTPVVRPAHISAFLQDRPIPEWCGVQHIHLSPSHHSGSKAASL' +
            'HWTSERVVSVLLLGLLPAAYLNPCSAMDYSLAAALTFMVTGALDKLLLTMFMGMPCRKLPRQG' +
            'FWHFQL*'),
            ('MLWRLSAVCGALGGRALLLRTPVVRPAHISAFLQDRPIPEWCGVQHIHLSPSHHSGSKAASL' +
            'HWTSERVVSVLLLGLLPAAYLNPCSAMYYSLAAALTLHGHWGLGQVVTDYVHGDALQKAAKAG' +
            'LLALSALTFAGLCYFNYHDVGICKPLPCCGSS*'),
            ('MAVLWRLSAVCGAPTARDRRPSSVASNSSGQTCSYLSISSGPTYPRMVWSAAHTLVTEPPFW' +
            'LQGCISPLD*')
        ]

        self.results = []
        for sample in self.sample:
            self.results.append(describe_protein(self.reference, sample))

    def _single_variant(self, sample, expected):
        """
        General single variant test.
        """
        description = describe.describe_protein(self.reference, sample)
        assert description[0].type == expected[0]
        assert description[0].start == expected[1]
        assert description[0].end == expected[2]
        assert description[0].sample_start == expected[3]
        assert description[0].sample_end == expected[4]
        assert description[0].deleted[0].sequence == expected[5]
        assert description[0].inserted[0].sequence == expected[6]
        assert str(description[0]) == expected[7]


    def test1(self):
        """
        Test 1.
        """
        assert str(self.results[0][0]) == '[Ser8_Ala9del;Leu101Phefs*34]'


    def test2(self):
        """
        Test 2.
        """
        assert (str(self.results[1][0]) ==
            '[Ala2_Val3del;Asp92Tyr;Ala152Profs*9]')


    def test3(self):
        """
        Test 3.
        """
        assert str(self.results[2][0]) == ('Leu14_Leu159delinsProThrAlaArg' +
            'AspArgArgProSerSerValAlaSerAsnSerSerGlyGlnThrCysSerTyrLeuSerIle' +
            'SerSerGlyProThrTyrProArgMetValTrpSerAlaAlaHisThrLeuValThrGluPro' +
            'ProPheTrpLeuGlnGlyCysIleSerProLeuAsp')


    def test4(self):
        """
        Test 4.
        """
        assert (self.results[0][0].nhgvs() ==
            '[8_9del;101_158delinsFMVTGALDKLLLTMFMGMPCRKLPRQGFWHFQ]')


    def test5(self):
        """
        Test 5.
        """
        assert (self.results[1][0].nhgvs() ==
            '[2_3del;92D>Y;152_159delinsPLPCCGSS]')


    def test6(self):
        """
        Test 6.
        """
        assert self.results[2][0].nhgvs() == ('14_159delinsPTARDRRPSSVASNSSG' +
            'QTCSYLSISSGPTYPRMVWSAAHTLVTEPPFWLQGCISPLD')


    def test7(self):
        """
        Test 7.
        """
        assert str(self.results[0][1]) == '101_133fs+2'


    def test8(self):
        """
        Test 8.
        """
        assert str(self.results[1][1]) == '152_159fs+2'


    def test9(self):
        """
        Test 9.
        """
        assert str(self.results[2][1]) == '14_67fs+1'


    def test10(self):
        """
        Test 10.
        """
        description = describe_protein('MDYSLAAALTLHGH',
            'MTIPWRSPHFHGH')
        assert (str(description) ==
            'Asp2_Leu11delinsThrIleProTrpArgSerProHisPhe')


#    def test5(self):
#        """
#        Test 5.
#        """
#        self._single_variant('ACGTCGATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
#            ('none', 0, 0, 0, 0, '', '', '='))
#
#
#    def test6(self):
#        """
#        Test 6.
#        """
#        self._single_variant('ACGTCGGTTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
#            ('subst', 7, 7, 7, 7, 'A', 'G', '7A>G'))
#
#
#    def test7(self):
#        """
#        Test 7.
#        """
#        self._single_variant('ACGTCGTTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
#            ('del', 7, 7, 6, 7, 'A', '', '7del'))
#
#
#    def test8(self):
#        """
#        Test 8.
#        """
#        self._single_variant('ACGTCGTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
#            ('del', 7, 8, 6, 7, 'AT', '', '7_8del'))
#
#
#    def test9(self):
#        """
#        Test 9.
#        """
#        self._single_variant('ACGTCGCATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
#            ('ins', 6, 7, 7, 7, '', 'C', '6_7insC'))
#
#
#    def test10(self):
#        """
#        Test 10.
#        """
#        self._single_variant('ACGTCGCCATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
#            ('ins', 6, 7, 7, 8, '', 'CC', '6_7insCC'))
#
#
#    def test11(self):
#        """
#        Test 11.
#        """
#        self._single_variant('ACGTCGAATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
#            ('dup', 7, 7, 8, 8, '', 'A', '7dup'))
#
#
#    def test12(self):
#        """
#        Test 12.
#        """
#        self._single_variant('ACGTCGAGATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
#            ('dup', 6, 7, 8, 9, '', 'GA', '6_7dup'))
#
#
#    def test13(self):
#        """
#        Test 13.
#        """
#        self._single_variant('ACGTCGACGATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
#            ('dup', 5, 7, 8, 10, '', 'CGA', '5_7dup'))
#
#
#    def test14(self):
#        """
#        Test 14.
#        """
#        self._single_variant('ACGTCGCGAATCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
#            ('inv', 7, 11, 7, 11, 'ATTCG', 'CGAAT', '7_11inv'))
#
#
#    def test15(self):
#        """
#        Test 15.
#        """
#        self._single_variant('ACGTCGCCTTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
#            ('delins', 7, 7, 7, 8, 'A', 'CC', '7delinsCC'))
#
#
#    def test16(self):
#        """
#        Test 16.
#        """
#        self._single_variant('ACGTCGATTCGCTAGCTTCGTTTTGATAGATAGAGATATAGAGAT',
#            ('delins', 21, 23, 21, 24, 'GGG', 'TTTT', '21_23delinsTTTT'))
#
#
#    def test17(self):
#        """
#        Test 17.
#        """
#        self._single_variant('ACGTCTCTTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
#            ('inv', 6, 7, 6, 7, 'GA', 'TC', '6_7inv'))
#
#
#    def test18(self):
#        """
#        Test 18.
#        """
#        self._single_variant('ACGTCGTCTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT',
#            ('delins', 7, 8, 7, 8, 'AT', 'TC', '7_8delinsTC'))
