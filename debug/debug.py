# debuging that loss is zero for one of the encoder ouput the derivative is also zero
params.loss_excluded_from_training
list(outputs.parameters())
for p in outputs.named_parameters():
    print(p)

cc = outputs.parameters().next()
