from create_matching_graph import Matching_graph

# from sim3d import Sim3D
# from qecsim.src.pauli import Pauli
from visualize import visualize
from layout import Hexagonal_layout

# from decoding_process import Decoding_process


g5 = Matching_graph(5)
d5_setup = Hexagonal_layout(5)
d7_setup = Hexagonal_layout(7)


def test_visualize_matching_graph():
    vis = visualize(3, d5_setup)
    vis.plot_matching_graph(matching_graph=g5)


# test_visualize_matching_graph()


def test_plot_syndrome():
    sim3d = Sim3D(5)
    decoder = Decoding_process(sim3d)
    error_list_X, error_list_Z = sim3d.create_error_list(0.01)
    error, synd_x_hist, synd_z_hist, error_history = sim3d.run_single_cycle_track(
        error_list_X, error_list_Z
    )

    vis = visualize(5, d5_setup)
    vis.plot_syndrome(synd_x_hist, synd_z_hist, error_history=error_history)


# test_plot_syndrome()


def test_plot_2d_layer():
    #    sim3d = Sim3D(5)
    # decoder = Decoding_process(sim3d)
    g5 = Matching_graph(5, n_layers=1)
    vis = visualize(5, d5_setup)
    # , 'red'))
    vis.plot_2d_layer(("green", "red"), plot_edges=True, matching_graph=g5)
    # plot_edges=True, matching_graph=g5)


test_plot_2d_layer()
