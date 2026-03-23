# imports
import numpy as np
from tangermeme.predict import predict
from tangermeme.ism import saturation_mutagenesis


def greedy_ism(
    model,
    X,
    seq_len=2114,
    ism_window=500,
    k=1,
    max_iter=10,
    verbose=False,
    mode="maximize",  # New parameter
):
    """
    Perform greedy ISM on a sequence.

    Parameters
    ----------
    model : BPNet
        The model to use for prediction.
    X : torch.Tensor
        The input sequence as a one-hot encoded tensor.
    seq_len : int
        The length of the input sequence.
    ism_window : int
        The window size for ISM.
    k : int
        The number of top positions to select.
    max_iter : int
        The maximum number of iterations to perform.
    verbose : bool
        Whether to print progress messages.
    mode : str
        Whether to 'maximize' or 'minimize' the model score. Default is 'maximize'.

    Returns
    -------
    X_ : torch.Tensor
        The modified sequence after greedy ISM.
    """

    # Define seq limits
    seq_center = seq_len // 2
    seq_start = seq_center - ism_window // 2
    seq_end = seq_center + ism_window // 2

    # Initial copy of the tensor
    X_ = X.clone()

    # Get the original prediction
    y0 = predict(model, X_).cpu().numpy()[0]
    print(f"Iteration 0/{max_iter} - y_hat: {y0[0]:.3f}")

    # Run ISM across rounds
    for i in range(max_iter):
        if verbose:
            print(f"Iteration {i+1}/{max_iter}", end="")

        # Run ISM
        y0, y_hat = saturation_mutagenesis(
            model,
            X_,
            start=seq_start,
            end=seq_end,
            raw_outputs=True,
            verbose=False,
        )

        # Get the updates
        attr = y_hat[:, :, :, 0] - y0[:, None, None, 0]
        attr = attr.numpy()

        # Take the average across the batch
        if X.shape[0] == 2:
            attr_mean = np.mean(attr, axis=0)
        else:
            attr_mean = attr.squeeze(0)

        # Flip attribution if minimizing
        if mode == "minimize":
            attr_mean = -attr_mean

        # Get the top k positions
        idx = np.argsort(attr_mean.ravel())[:-k-1:-1]
        indices = np.column_stack(np.unravel_index(idx, attr_mean.shape))[0]

        # if verbose print out the new y_hat value
        if verbose:
            y_hat_np = y_hat[:, :, :, 0].cpu().numpy()
            print(f" - y_hat: {y_hat_np[0, indices[0], indices[1]]:.3f}")

        # Get the edit channel and position
        edit_channel = indices[0]
        edit_position = seq_start + indices[1]

        # zero out the position across the channels
        X_[0, :, edit_position] = 0
        X_[0, edit_channel, edit_position] = 1

    return X_
