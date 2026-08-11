import unittest
from main import GF256, GF256Poly, RSCodec

class TestReedSolomonDecode(unittest.TestCase):
    def setUp(self):
        self.codec = RSCodec(255,251)
        self.msg = "this is a test message"
        self.codeword = self.codec.encode(self.msg)

    def test_zero_syndromes_for_uncorrupted_message(self):
        synds = self.codec.decode(self.codeword.coeffs)

        self.assertEqual(
            synds,
            [0] * self.codec.nsym,
            f"Expected all zeros, got {synds}",
        )

    def test_more_syndromes_for_corrupted_message(self):

        encoded_msg = bytes(self.codeword.coeffs)

        corrupted_msg = bytearray(encoded_msg)
        corrupted_msg[0] ^= 0xFF

        synds = self.codec.decode(corrupted_msg)
        n = sum(synds)

        self.assertTrue(
            any(s != 0 for s in synds),
            "Syndromes should be non-zero when errors are present",
        )

    def test_decode_matches_via_horners_method(self):
        encoded_msg = bytes(self.codeword.coeffs)

if __name__ == "__main__":
    unittest.main()
