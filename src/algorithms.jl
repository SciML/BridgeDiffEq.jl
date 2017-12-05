abstract type BridgeAlgorithm <: DEAlgorithm end
struct BridgeEuler <: BridgeAlgorithm end
struct BridgeHeun <: BridgeAlgorithm end
struct BridgeSRK <: BridgeAlgorithm end
struct BridgeR3 <: BridgeAlgorithm end
struct BridgeBS3 <: BridgeAlgorithm end
