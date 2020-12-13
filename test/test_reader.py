from snpfrag.core.fsareader import FsaReader


def test_reader():

    reader = FsaReader("examples/A01.fsa")

    assert reader.dye_names == ['6-FAM', 'HEX', 'NED', 'ROX']
    assert "DATA1" in reader.channel_names
    size_channel = reader.channel_from_dye_name("ROX") 
    assert size_channel == "DATA4"

    ROX_SIZE = [35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490 , 500 ]

    reader.compute_regression(size_channel,ROX_SIZE)
    assert reader.regression.stderr < 0.005



    
